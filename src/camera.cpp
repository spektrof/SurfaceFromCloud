#include "camera.h"
#include <math.h>
#include <algorithm>
#include <QtMath>

// Initializes a new instance of the class.
Camera::Camera(void) : m_eye(0.0f, -20.0f, 20.0f), m_at(0.0f, 0.0f, 0.0f), m_up(0.0f, 1.0f, 0.0f), m_speed(16.0f), m_goFw(0), m_goRight(0), m_slow(false)
{
	ortho_rect = QVector4D(-15.0f, 15.0f, 15.0f, -15.0f);
	nearPlane = -505.0f;
	farPlane = 8000.0f;

	dir = 1;
	dir_up = 1;
	zoom = 1.0f;

	radius = 10.0f;
	theta = 0.5f;
	phi = 0.5f;
	m_eye = QVector3D(radius * cos(theta * PI) * sin(phi * PI), radius * cos(phi * PI), radius * sin(theta * PI) * sin(phi * PI));

	m_matProj.setToIdentity();
	m_matProj.ortho(ortho_rect.x(), ortho_rect.y(), ortho_rect.z(), ortho_rect.w(), nearPlane, farPlane);
	m_viewMatrix.setToIdentity();
	m_viewMatrix.lookAt(m_eye, m_at, m_up);
	m_matViewProj = m_matProj * m_viewMatrix;
	//SetView(m_eye,m_at,m_up);

	radius_up = sqrt(1 + radius * radius);
	phi_up = 0.5f - asinf(1.0f / radius) / PI;
	theta_up = 0.5f;
	m_up = QVector3D(radius_up * cos(theta_up * PI) * sin(phi_up * PI), radius_up * cos(phi_up * PI), radius_up * sin(theta_up * PI) * sin(phi_up * PI));

	QVector3D tmp = m_eye;
	tmp.normalize();
	m_up = m_up - (radius * tmp);
	m_up.normalize();

//	qDebug() << "UP VECTOR: " << m_up << "\n";

}

Camera::Camera(QVector3D _eye, QVector3D _at, QVector3D _up) : m_speed(16.0f), m_goFw(0), m_goRight(0), m_dist(10), m_slow(false)
{
	SetView(_eye, _at, _up);
}

Camera::~Camera(void)
{
}

void Camera::SetView(QVector3D _eye, QVector3D _at, QVector3D _up)
{
	//m_eye = _eye;
	//m_at = _at;
	//m_up = _up;

	/*m_fw = m_at - m_eye;
	m_fw.normalize();

	m_st = (QVector3D::crossProduct(m_fw, m_up));
	m_st.normalize();

	m_dist = (m_at - m_eye).length();

	m_u = atan2f(m_fw.z(), m_fw.x());
	m_v = acosf(m_fw.y());*/
}

void Camera::SetProj(float _angle, float _aspect, float _zn, float _zf, float a, float b)
{
	m_matProj.setToIdentity();
	m_matProj.ortho(_angle, _aspect, _zn, _zf, a, b);
	m_matViewProj = m_matProj * m_viewMatrix;
}

QMatrix4x4 Camera::GetViewMatrix()
{
	return m_viewMatrix;
}

bool Camera::Update(float _deltaTime)
{
	if (m_goFw == 0 && m_goRight == 0) return false;
	
	phi += dir * (m_goFw* 0.05)*m_speed*_deltaTime;
	phi_up += dir_up * (m_goFw* 0.05)*m_speed*_deltaTime;
	theta += (m_goRight* 0.05)*m_speed*_deltaTime;
	theta_up += (m_goRight* 0.05)*m_speed*_deltaTime;

	if ( m_goRight && (theta > 2.0f || theta < 0))	//ezek okoznak még hibát..., pontosabban a 2.0f - emiatt LOG BENN MARAD
	{
		theta = theta > 2.0f ? 0.0f : 2.0f;	
	}

	if (m_goRight && (theta_up > 2.0f || theta_up < 0))
	{
		theta_up = theta_up > 2.0f ? 0.0f : 2.0f;
	}

	if (phi > 1.0f || phi <0)
	{
		dir *= -1;
		phi = (int)phi;
		if (theta >= 1.0f)
		{
			if (theta >= 1.5f)
			{
				theta = 1.0f - (2.0f - theta);
			}
			else
				theta = theta - 1.0f;
		}
		else
		{
			if (theta >= 0.5f)
			{
				theta = 2.0f - (1.0f - theta);
			}
			else
				theta = 1.0f + theta_up;
		}
		
	}

	if (phi_up > 1.0f || phi_up < 0)
	{
		//theta_up = 1.0f + (theta_up > 1 ? -1 * (theta_up - 1.0f) : 1 * (1.0f - theta_up));
		dir_up *= -1;
		phi_up = (int)phi_up;
		if (theta_up >= 1.0f)
		{
			if (theta_up >= 1.5f)
			{
				theta_up = 1.0f - (2.0f - theta_up);	//jó
			}
			else
				theta_up = theta_up - 1.0f;
		}
		else
		{
			if (theta_up >= 0.5f)
			{
				theta_up = 2.0f - (1.0f - theta_up);	//jó
			}
			else
				theta_up = 1.0f + theta_up;
		}
	}

	qDebug() << phi << "\n";
	qDebug() << theta << "\n";

	m_up = QVector3D(radius_up * cos(theta_up * PI) * sin(phi_up * PI), radius_up * cos(phi_up * PI), radius_up * sin(theta_up * PI) * sin(phi_up * PI));
	m_eye = QVector3D(radius * cos(theta * PI) * sin(phi * PI), radius * cos(phi * PI), radius * sin(theta * PI) * sin(phi * PI));

	QVector3D tmp = m_eye;
	tmp.normalize();
	m_up = m_up - (radius * tmp);
	m_up.normalize();

	qDebug() << m_eye << "\n";
	qDebug() << m_up << "\n";

	m_viewMatrix.setToIdentity();
	m_viewMatrix.lookAt(m_eye, m_at, m_up);
	m_matViewProj = m_matProj * m_viewMatrix;
	return true;
}

template<typename T>
T clamp(T x, T minVal, T maxVal)
{
	return std::min(std::max(x, minVal), maxVal);
}

void Camera::UpdateUV(float du, float dv)
{
	m_u += du;
	m_v = clamp<float>( m_v + dv, 0.1f, 3.1f);

	m_at = m_eye + m_dist * QVector3D(cosf(m_u)*sinf(m_v),
		cosf(m_v),
		sinf(m_u)*sinf(m_v));

	m_fw = m_at - m_eye;
	m_fw.normalize();

	m_st = QVector3D::crossProduct(m_fw, m_up);
	m_st.normalize();
}

void Camera::SetSpeed(float _val)
{
	m_speed = _val;
}

void Camera::Resize(int _w, int _h)
{
	
}

void Camera::KeyboardDown(QKeyEvent *event)
{
	switch (event->key())
	{
	case Qt::Key_W:
	{
		m_goFw = 1;
		break;
	}
	case Qt::Key_S:
	{
		m_goFw = -1;
		break;
	}
	case Qt::Key_A:
	{
		m_goRight = -1;
		break;
	}
	case Qt::Key_D:
	{
		m_goRight = 1;
		break;
	}
	/*case SDLK_LSHIFT:
	case SDLK_RSHIFT:
		if (!m_slow)
		{
			m_slow = true;
			m_speed /= 4.0f;
		}
		break;
	case SDLK_w:
		m_goFw = 1;
		break;
	case SDLK_s:
		m_goFw = -1;
		break;
	case SDLK_a:
		m_goRight = -1;
		break;
	case SDLK_d:
		m_goRight = 1;
		break;*/
	}
}

void Camera::KeyboardUp(QKeyEvent *event)
{
	float current_speed = m_speed;
	switch (event->key())
	{
	case Qt::Key_W:
	{
		m_goFw = 0;
		break;
	}
	case Qt::Key_S:
	{
		m_goFw = 0;
		break;
	}
	case Qt::Key_A:
	{

		m_goRight = 0;
		break;
	}
	case Qt::Key_D:
	{

		m_goRight = 0;
		break;
	}
	/*case SDLK_LSHIFT:
	case SDLK_RSHIFT:
		if (m_slow)
		{
			m_slow = false;
			m_speed *= 4.0f;
		}
		break;
	case SDLK_w:
	case SDLK_s:
		m_goFw = 0;
		break;
	case SDLK_a:
	case SDLK_d:
		m_goRight = 0;
		break;*/
	}
}

void Camera::wheelEvent(QWheelEvent *event)
{
	float magnitude = event->angleDelta().y() / 240.0f;

	zoom = zoom + magnitude > 0 ? zoom + magnitude : zoom;
	SetProj( ortho_rect.x() / zoom, ortho_rect.y() / zoom, ortho_rect.z() / zoom, ortho_rect.w() / zoom,  nearPlane, farPlane);
}

/*
void Camera::MouseMove(SDL_MouseMotionEvent& mouse)
{
	if (mouse.state & SDL_BUTTON_LMASK)
	{
		UpdateUV(mouse.xrel / 100.0f, mouse.yrel / 100.0f);
	}
}
*/
void Camera::LookAt(QVector3D _at)
{
	SetView(m_eye, _at, m_up);
}

