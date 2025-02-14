OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1887551) q[0];
sx q[0];
rz(-2.3030757) q[0];
sx q[0];
rz(1.4128348) q[0];
rz(0.38263327) q[1];
sx q[1];
rz(1.1525947) q[1];
sx q[1];
rz(11.424185) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.45158) q[0];
sx q[0];
rz(-0.40209189) q[0];
sx q[0];
rz(1.3908338) q[0];
x q[1];
rz(0.76337645) q[2];
sx q[2];
rz(-0.65033153) q[2];
sx q[2];
rz(-1.842832) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9153629) q[1];
sx q[1];
rz(-2.1945476) q[1];
sx q[1];
rz(-0.28304328) q[1];
rz(-pi) q[2];
rz(-0.88291165) q[3];
sx q[3];
rz(-1.7241244) q[3];
sx q[3];
rz(-0.49082498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4932058) q[2];
sx q[2];
rz(-1.1492665) q[2];
sx q[2];
rz(-3.1374078) q[2];
rz(-0.74797136) q[3];
sx q[3];
rz(-2.3145521) q[3];
sx q[3];
rz(0.23676693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0325322) q[0];
sx q[0];
rz(-0.91962686) q[0];
sx q[0];
rz(0.064706651) q[0];
rz(-2.0789355) q[1];
sx q[1];
rz(-0.92667842) q[1];
sx q[1];
rz(-2.0978755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4499265) q[0];
sx q[0];
rz(-1.1083034) q[0];
sx q[0];
rz(-1.3426379) q[0];
x q[1];
rz(1.6690786) q[2];
sx q[2];
rz(-1.543902) q[2];
sx q[2];
rz(-1.310845) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5001561) q[1];
sx q[1];
rz(-1.078633) q[1];
sx q[1];
rz(-1.3970966) q[1];
rz(-pi) q[2];
rz(-1.4426368) q[3];
sx q[3];
rz(-1.7944873) q[3];
sx q[3];
rz(2.8552289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8740497) q[2];
sx q[2];
rz(-2.0127313) q[2];
sx q[2];
rz(-1.4675325) q[2];
rz(-2.2360146) q[3];
sx q[3];
rz(-2.0445243) q[3];
sx q[3];
rz(-2.9300698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8242789) q[0];
sx q[0];
rz(-2.2792094) q[0];
sx q[0];
rz(-0.30340075) q[0];
rz(-0.41269451) q[1];
sx q[1];
rz(-1.7957325) q[1];
sx q[1];
rz(2.5491098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2148662) q[0];
sx q[0];
rz(-2.0030372) q[0];
sx q[0];
rz(-1.809948) q[0];
rz(-pi) q[1];
rz(0.5634339) q[2];
sx q[2];
rz(-1.7220925) q[2];
sx q[2];
rz(-0.39416692) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35960782) q[1];
sx q[1];
rz(-2.4949412) q[1];
sx q[1];
rz(2.9293865) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72393317) q[3];
sx q[3];
rz(-1.9125126) q[3];
sx q[3];
rz(-0.011499876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5182284) q[2];
sx q[2];
rz(-0.48064226) q[2];
sx q[2];
rz(2.9474958) q[2];
rz(0.57191166) q[3];
sx q[3];
rz(-1.5827936) q[3];
sx q[3];
rz(1.5552049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4647575) q[0];
sx q[0];
rz(-2.0763626) q[0];
sx q[0];
rz(1.7179426) q[0];
rz(-2.8061197) q[1];
sx q[1];
rz(-0.60233855) q[1];
sx q[1];
rz(2.4981892) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62221974) q[0];
sx q[0];
rz(-1.046139) q[0];
sx q[0];
rz(-0.93726802) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57482608) q[2];
sx q[2];
rz(-1.4919623) q[2];
sx q[2];
rz(-3.0835033) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7285541) q[1];
sx q[1];
rz(-0.87552363) q[1];
sx q[1];
rz(2.1785582) q[1];
rz(-1.8845064) q[3];
sx q[3];
rz(-2.6711552) q[3];
sx q[3];
rz(1.1251118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8947577) q[2];
sx q[2];
rz(-2.0474696) q[2];
sx q[2];
rz(2.3168054) q[2];
rz(0.60683933) q[3];
sx q[3];
rz(-1.9347128) q[3];
sx q[3];
rz(1.8230033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7428699) q[0];
sx q[0];
rz(-2.5720808) q[0];
sx q[0];
rz(-2.8745765) q[0];
rz(2.6737402) q[1];
sx q[1];
rz(-2.0303969) q[1];
sx q[1];
rz(1.9652479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0396966) q[0];
sx q[0];
rz(-1.5202664) q[0];
sx q[0];
rz(0.021357892) q[0];
rz(-pi) q[1];
rz(-0.68525542) q[2];
sx q[2];
rz(-0.90327679) q[2];
sx q[2];
rz(-1.670862) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49217192) q[1];
sx q[1];
rz(-2.6876059) q[1];
sx q[1];
rz(-1.6976056) q[1];
x q[2];
rz(0.72132206) q[3];
sx q[3];
rz(-2.0287987) q[3];
sx q[3];
rz(1.3617226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0331369) q[2];
sx q[2];
rz(-1.5507853) q[2];
sx q[2];
rz(2.1948658) q[2];
rz(-1.6999792) q[3];
sx q[3];
rz(-1.0333034) q[3];
sx q[3];
rz(-1.3942963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93447584) q[0];
sx q[0];
rz(-1.2756791) q[0];
sx q[0];
rz(-1.7913272) q[0];
rz(-2.4840202) q[1];
sx q[1];
rz(-1.5788199) q[1];
sx q[1];
rz(-1.2616166) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5367071) q[0];
sx q[0];
rz(-1.1632445) q[0];
sx q[0];
rz(-0.85257963) q[0];
rz(-2.4305651) q[2];
sx q[2];
rz(-2.3287352) q[2];
sx q[2];
rz(2.6136398) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.99743479) q[1];
sx q[1];
rz(-0.23845928) q[1];
sx q[1];
rz(2.8644117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2618939) q[3];
sx q[3];
rz(-2.8692103) q[3];
sx q[3];
rz(2.4177944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6235846) q[2];
sx q[2];
rz(-2.035718) q[2];
sx q[2];
rz(-0.27274954) q[2];
rz(2.2128211) q[3];
sx q[3];
rz(-2.2835968) q[3];
sx q[3];
rz(-1.6509008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6269094) q[0];
sx q[0];
rz(-1.0310443) q[0];
sx q[0];
rz(2.2169901) q[0];
rz(2.5968016) q[1];
sx q[1];
rz(-1.3489312) q[1];
sx q[1];
rz(-3.0371688) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69014617) q[0];
sx q[0];
rz(-1.4826097) q[0];
sx q[0];
rz(-3.0549391) q[0];
rz(1.7288389) q[2];
sx q[2];
rz(-1.9812532) q[2];
sx q[2];
rz(2.4685658) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.663781) q[1];
sx q[1];
rz(-2.7594477) q[1];
sx q[1];
rz(-2.2078329) q[1];
x q[2];
rz(-0.84568328) q[3];
sx q[3];
rz(-0.50456797) q[3];
sx q[3];
rz(-2.9364236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8203848) q[2];
sx q[2];
rz(-1.2106004) q[2];
sx q[2];
rz(1.4804117) q[2];
rz(-0.004247578) q[3];
sx q[3];
rz(-0.85631266) q[3];
sx q[3];
rz(-0.64928865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55411196) q[0];
sx q[0];
rz(-0.09621796) q[0];
sx q[0];
rz(1.2531248) q[0];
rz(-0.15086497) q[1];
sx q[1];
rz(-2.3020709) q[1];
sx q[1];
rz(-1.2996659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4010581) q[0];
sx q[0];
rz(-2.5098233) q[0];
sx q[0];
rz(2.6982672) q[0];
x q[1];
rz(2.9787709) q[2];
sx q[2];
rz(-2.0042684) q[2];
sx q[2];
rz(-2.2107645) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35132699) q[1];
sx q[1];
rz(-2.4557487) q[1];
sx q[1];
rz(0.02355656) q[1];
rz(-pi) q[2];
rz(1.9261041) q[3];
sx q[3];
rz(-1.7607712) q[3];
sx q[3];
rz(3.0440898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3018939) q[2];
sx q[2];
rz(-1.3816741) q[2];
sx q[2];
rz(1.985644) q[2];
rz(0.65129748) q[3];
sx q[3];
rz(-0.39172253) q[3];
sx q[3];
rz(-2.7294066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10525178) q[0];
sx q[0];
rz(-1.3407433) q[0];
sx q[0];
rz(-2.3403008) q[0];
rz(-2.2835412) q[1];
sx q[1];
rz(-0.67887226) q[1];
sx q[1];
rz(-0.37989315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86738619) q[0];
sx q[0];
rz(-2.634122) q[0];
sx q[0];
rz(1.9486289) q[0];
rz(-0.97486603) q[2];
sx q[2];
rz(-2.9306378) q[2];
sx q[2];
rz(0.67529222) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8752138) q[1];
sx q[1];
rz(-2.3734178) q[1];
sx q[1];
rz(-1.1298864) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5809459) q[3];
sx q[3];
rz(-2.0119609) q[3];
sx q[3];
rz(-0.084189296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6259367) q[2];
sx q[2];
rz(-2.8069324) q[2];
sx q[2];
rz(-0.77896172) q[2];
rz(-0.37911478) q[3];
sx q[3];
rz(-1.4232676) q[3];
sx q[3];
rz(0.085722119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871386) q[0];
sx q[0];
rz(-1.7310646) q[0];
sx q[0];
rz(2.6692303) q[0];
rz(0.27913276) q[1];
sx q[1];
rz(-1.0529073) q[1];
sx q[1];
rz(-0.36769029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5147509) q[0];
sx q[0];
rz(-1.923133) q[0];
sx q[0];
rz(0.3120089) q[0];
x q[1];
rz(-0.39989465) q[2];
sx q[2];
rz(-1.7755819) q[2];
sx q[2];
rz(2.9861394) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93802261) q[1];
sx q[1];
rz(-1.9782776) q[1];
sx q[1];
rz(1.8304871) q[1];
rz(-pi) q[2];
rz(0.69766694) q[3];
sx q[3];
rz(-0.945745) q[3];
sx q[3];
rz(2.5755957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3231861) q[2];
sx q[2];
rz(-2.0995993) q[2];
sx q[2];
rz(2.6673178) q[2];
rz(-0.91931528) q[3];
sx q[3];
rz(-1.5761458) q[3];
sx q[3];
rz(-1.3924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267088) q[0];
sx q[0];
rz(-1.8395431) q[0];
sx q[0];
rz(-1.607847) q[0];
rz(2.907091) q[1];
sx q[1];
rz(-2.0461743) q[1];
sx q[1];
rz(2.8428427) q[1];
rz(-1.7099829) q[2];
sx q[2];
rz(-2.5659106) q[2];
sx q[2];
rz(2.8775712) q[2];
rz(-1.9145253) q[3];
sx q[3];
rz(-3.0889441) q[3];
sx q[3];
rz(-2.3788388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
