#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QVector3D>
#include <QMatrix4x4>
#include <QKeyEvent>


class Camera
{
public:
	Camera(void);
	Camera(QVector3D _eye, QVector3D _at, QVector3D _up);
	~Camera(void);

	QMatrix4x4 GetViewMatrix();

	bool Update(float _deltaTime);

	void SetView(QVector3D _eye, QVector3D _at, QVector3D _up);
	void SetProj(float _angle, float _aspect, float _zn, float _zf, float a = 0.0f, float b = 0.0f);
	void LookAt(QVector3D _at);

	void SetSpeed(float _val);
	QVector3D GetEye() { return m_eye; }

	QVector3D GetAt() { return m_at; }

	QVector3D GetUp() { return m_up; }

	QMatrix4x4 GetProj() { return m_matProj; }

	QMatrix4x4 GetViewProj() {	return m_matViewProj; }

	void Resize(int _w, int _h);

	void KeyboardDown(QKeyEvent *event);
	void KeyboardUp(QKeyEvent *event);
	void wheelEvent(QWheelEvent *event);
	/*void MouseMove(SDL_MouseMotionEvent& mouse);*/

private:
#define PI 3.14159265359

	void UpdateUV(float du, float dv);

	bool	m_slow;

	// The camera position.
 
	// The u spherical coordinate of the spherical coordinate pair (u,v) denoting the
	// current viewing direction from the view position m_eye. 
 
	float	m_u;

	// <summary>
	// The v spherical coordinate of the spherical coordinate pair (u,v) denoting the
	// current viewing direction from the view position m_eye. 
 
	float	m_v;

	// The distance of the look at point from the camera. 
	float	m_dist;

	// The unit vector pointing towards the viewing direction.
	QVector3D	m_fw;
	// The unit vector pointing to the 'right'
	QVector3D	m_st;

	QMatrix4x4	m_matProj;

	float	m_goFw;
	float	m_goRight;

	//---------------------
	QVector4D ortho_rect;

	float nearPlane, farPlane;
	float radius, radius_up;
	float theta, phi, phi_up, theta_up;
	int dir, dir_up;
	float zoom;

	QVector3D	m_eye;
	QVector3D	m_up;
	QVector3D	m_at;

	float		m_speed;

	QMatrix4x4	m_viewMatrix;
	QMatrix4x4	m_matViewProj;

};

