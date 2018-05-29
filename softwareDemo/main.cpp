#include <GL/freeglut.h>
#include <GL/gl.h>
#include <GL/glu.h>

struct GLpoint
{
    GLfloat x,y;
};

inline void drawDot(GLint x,GLint y)
{
    glBegin(GL_POINTS);
        glVertex2i(x,y);
    glEnd();
}

void myInit()
{
    glClearColor(1.0,1.0,1.0,0.0);//设置背景颜色为亮白
    glColor3f(0.0f,0.0f,0.0f);//设置绘图颜色为黑色
    glPointSize(4.0);//设置点的大小为4x4像素
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0,640.0,0.0,480.0);
}

void myDisplay()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_POINTS);
        glVertex2i(100,50);
        glVertex2i(100,130);
        glVertex2i(150,130);
    glEnd();
    glFlush();
}

int main(int argc,char **argv)
{
    glutInit(&argc,argv);//初始化工具包
    glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);//设置显示模式
    glutInitWindowSize(640,480);//设置窗口大小
    glutInitWindowPosition(100,100);//设置屏幕上窗口位置
    glutCreateWindow("my first attemot");//打开带标题的窗口
    glutDisplayFunc(myDisplay);//注册重画回调函数
    myInit();
    glutMainLoop();
    return 0;
}