#version 130

in vec3 vs_out_pos;
out vec4 fs_out_col;

uniform mat4 viewProj;
uniform mat4 viewIprojI;
uniform mat4 modelI;
uniform mat4 model;
uniform mat4 view;

//uniform float r = 0.05f;

uniform int size = 0;
uniform vec3 centers[100];
uniform float radius[100];

uniform vec3 color = glm::vec3(0,0,1);
//uniform vec3 colors[100];	//later mybe

void getRay(in vec3 inVec, out vec3 rayOrig, out vec3 rayDir)
{
	// a világKR-ben a kattintásnak a közeli vágósíkon megfeleltetett pont koordinátái
	vec4 nearPt = viewIprojI * vec4(inVec.xy,-1, 1);
	// a világKR-ben a kattintásnak a távoli vágósíkon megfeleltetett pont koordinátái
	vec4 farPt  = viewIprojI * vec4(inVec.xy, 1, 1);

	rayOrig = nearPt.xyz/nearPt.w;

	vec3 rayEnd = farPt.xyz/farPt.w;
	rayDir  = normalize( rayEnd - rayOrig  );
}

float sphere_distance(vec3 c, vec3 x, float r)
{
	vec3 d = c - x;
	return dot(d,d) - r*r;
}

void main()
{
	if (size == 0) discard;

	vec3 rayOrig, rayDir;

	getRay(vs_out_pos, rayOrig, rayDir);

	rayOrig = (modelI * vec4(rayOrig, 1) ).xyz;
	rayDir  = (modelI * vec4(rayDir,  0) ).xyz;

	const int STEP_COUNT = 1750;
	const float STEP_LENGTH = 0.015f;

	//int first = -1;
	//int last_step = 0;

	for (int i=0; i<STEP_COUNT; ++i)
	{
		vec3 p = rayOrig + i*STEP_LENGTH*rayDir;
		for (int j=0; j< size; ++j)
		{
			if (sphere_distance(centers[j],p, radius[j]) <= 0)
			{
				//if (first == -1 ) 
				//{
					fs_out_col = vec4( color, 1/sqrt(3));
				//	first = j;
				//}
				//last_step = i;
				return;
			}
		}
		
	}

	/*if (first >= 0)
	{
		vec3 intersectionPoint = rayOrig + last_step*STEP_LENGTH*rayDir;

		vec4 clipPos = viewProj * vec4( intersectionPoint, 1 );

		float zndc = clipPos.z / clipPos.w; 
	
		float n = gl_DepthRange.near;
		float f = gl_DepthRange.far;
	
		gl_FragDepth = (f-n)/2 * zndc + (f+n)/2;
		fs_out_col = vec4( gl_FragDepth, gl_FragDepth, gl_FragDepth, 1/sqrt(3));

	}
	else*/
	
	discard;

	//--------------------------- MYBE need the depth data later
	/*
	// ha mélységet is akarsz hozzárendelni:
	vec3 intersectionPoint = rayOrig + threshold_reach[i]*STEP_LENGTH*rayDir;

	vec4 clipPos = viewProj * vec4( intersectionPoint, 1 );

	float zndc = clipPos.z / clipPos.w; 
	
	float n = gl_DepthRange.near;
	float f = gl_DepthRange.far;
	
	gl_FragDepth = (f-n)/2 * zndc + (f+n)/2;
	*/
	
}