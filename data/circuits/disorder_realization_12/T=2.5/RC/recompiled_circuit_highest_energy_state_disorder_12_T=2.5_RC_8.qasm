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
rz(-1.9528376) q[0];
sx q[0];
rz(-0.83851695) q[0];
sx q[0];
rz(1.7287579) q[0];
rz(0.38263327) q[1];
sx q[1];
rz(-1.9889979) q[1];
sx q[1];
rz(-1.9994073) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6467616) q[0];
sx q[0];
rz(-1.1755623) q[0];
sx q[0];
rz(-3.0656205) q[0];
rz(-pi) q[1];
rz(-0.76337645) q[2];
sx q[2];
rz(-0.65033153) q[2];
sx q[2];
rz(1.842832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22622977) q[1];
sx q[1];
rz(-2.1945476) q[1];
sx q[1];
rz(0.28304328) q[1];
x q[2];
rz(0.19742404) q[3];
sx q[3];
rz(-2.2490778) q[3];
sx q[3];
rz(1.9367644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6483868) q[2];
sx q[2];
rz(-1.1492665) q[2];
sx q[2];
rz(-0.0041848103) q[2];
rz(-2.3936213) q[3];
sx q[3];
rz(-2.3145521) q[3];
sx q[3];
rz(2.9048257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10906049) q[0];
sx q[0];
rz(-2.2219658) q[0];
sx q[0];
rz(-0.064706651) q[0];
rz(1.0626571) q[1];
sx q[1];
rz(-0.92667842) q[1];
sx q[1];
rz(1.0437171) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2118156) q[0];
sx q[0];
rz(-2.6295595) q[0];
sx q[0];
rz(-0.42590745) q[0];
rz(-pi) q[1];
rz(1.30322) q[2];
sx q[2];
rz(-3.0397085) q[2];
sx q[2];
rz(-2.6153878) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5001561) q[1];
sx q[1];
rz(-2.0629597) q[1];
sx q[1];
rz(-1.3970966) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4426368) q[3];
sx q[3];
rz(-1.7944873) q[3];
sx q[3];
rz(0.28636378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26754293) q[2];
sx q[2];
rz(-1.1288613) q[2];
sx q[2];
rz(1.4675325) q[2];
rz(-0.90557805) q[3];
sx q[3];
rz(-2.0445243) q[3];
sx q[3];
rz(-0.21152285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31731376) q[0];
sx q[0];
rz(-0.86238328) q[0];
sx q[0];
rz(-0.30340075) q[0];
rz(-2.7288981) q[1];
sx q[1];
rz(-1.7957325) q[1];
sx q[1];
rz(-2.5491098) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45771407) q[0];
sx q[0];
rz(-1.787583) q[0];
sx q[0];
rz(-0.44332645) q[0];
rz(0.27806313) q[2];
sx q[2];
rz(-0.5812656) q[2];
sx q[2];
rz(2.1991625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.35960782) q[1];
sx q[1];
rz(-0.64665142) q[1];
sx q[1];
rz(-2.9293865) q[1];
rz(-2.0140225) q[3];
sx q[3];
rz(-0.89689287) q[3];
sx q[3];
rz(1.2942838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5182284) q[2];
sx q[2];
rz(-2.6609504) q[2];
sx q[2];
rz(2.9474958) q[2];
rz(2.569681) q[3];
sx q[3];
rz(-1.5587991) q[3];
sx q[3];
rz(1.5552049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67683515) q[0];
sx q[0];
rz(-1.06523) q[0];
sx q[0];
rz(-1.42365) q[0];
rz(2.8061197) q[1];
sx q[1];
rz(-0.60233855) q[1];
sx q[1];
rz(0.64340341) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59601027) q[0];
sx q[0];
rz(-1.032858) q[0];
sx q[0];
rz(0.62278231) q[0];
rz(2.5667666) q[2];
sx q[2];
rz(-1.4919623) q[2];
sx q[2];
rz(-3.0835033) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.57694084) q[1];
sx q[1];
rz(-1.1168861) q[1];
sx q[1];
rz(2.3481525) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8845064) q[3];
sx q[3];
rz(-2.6711552) q[3];
sx q[3];
rz(2.0164808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.246835) q[2];
sx q[2];
rz(-1.094123) q[2];
sx q[2];
rz(-2.3168054) q[2];
rz(-0.60683933) q[3];
sx q[3];
rz(-1.2068799) q[3];
sx q[3];
rz(-1.3185893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3987228) q[0];
sx q[0];
rz(-0.56951183) q[0];
sx q[0];
rz(2.8745765) q[0];
rz(2.6737402) q[1];
sx q[1];
rz(-2.0303969) q[1];
sx q[1];
rz(-1.1763447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0396966) q[0];
sx q[0];
rz(-1.6213263) q[0];
sx q[0];
rz(-3.1202348) q[0];
x q[1];
rz(-0.77645923) q[2];
sx q[2];
rz(-1.0506223) q[2];
sx q[2];
rz(-2.7732244) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.63309489) q[1];
sx q[1];
rz(-2.0208686) q[1];
sx q[1];
rz(-0.061636713) q[1];
rz(-2.4202706) q[3];
sx q[3];
rz(-1.112794) q[3];
sx q[3];
rz(-1.3617226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0331369) q[2];
sx q[2];
rz(-1.5908073) q[2];
sx q[2];
rz(0.94672686) q[2];
rz(1.6999792) q[3];
sx q[3];
rz(-1.0333034) q[3];
sx q[3];
rz(1.3942963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(0.93447584) q[0];
sx q[0];
rz(-1.2756791) q[0];
sx q[0];
rz(-1.7913272) q[0];
rz(0.65757242) q[1];
sx q[1];
rz(-1.5788199) q[1];
sx q[1];
rz(1.879976) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5091357) q[0];
sx q[0];
rz(-2.2194891) q[0];
sx q[0];
rz(0.52059569) q[0];
x q[1];
rz(0.71102755) q[2];
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
rz(0.71260604) q[1];
sx q[1];
rz(-1.7999876) q[1];
sx q[1];
rz(1.6372174) q[1];
rz(-pi) q[2];
rz(0.17619074) q[3];
sx q[3];
rz(-1.7796081) q[3];
sx q[3];
rz(-3.1274019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6235846) q[2];
sx q[2];
rz(-1.1058747) q[2];
sx q[2];
rz(-2.8688431) q[2];
rz(0.92877156) q[3];
sx q[3];
rz(-0.85799587) q[3];
sx q[3];
rz(1.4906918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(1.6729148) q[0];
sx q[0];
rz(-0.12355655) q[0];
sx q[0];
rz(0.79609032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4127537) q[2];
sx q[2];
rz(-1.9812532) q[2];
sx q[2];
rz(0.67302683) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1510091) q[1];
sx q[1];
rz(-1.2663453) q[1];
sx q[1];
rz(-0.23465381) q[1];
x q[2];
rz(-2.7905043) q[3];
sx q[3];
rz(-1.9410053) q[3];
sx q[3];
rz(2.5552487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8203848) q[2];
sx q[2];
rz(-1.2106004) q[2];
sx q[2];
rz(1.661181) q[2];
rz(3.1373451) q[3];
sx q[3];
rz(-0.85631266) q[3];
sx q[3];
rz(-0.64928865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5874807) q[0];
sx q[0];
rz(-0.09621796) q[0];
sx q[0];
rz(-1.8884678) q[0];
rz(-0.15086497) q[1];
sx q[1];
rz(-0.8395218) q[1];
sx q[1];
rz(-1.8419267) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74053451) q[0];
sx q[0];
rz(-0.63176934) q[0];
sx q[0];
rz(-0.44332544) q[0];
rz(-pi) q[1];
rz(-2.0093654) q[2];
sx q[2];
rz(-1.718443) q[2];
sx q[2];
rz(-2.4327337) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38176259) q[1];
sx q[1];
rz(-0.88517939) q[1];
sx q[1];
rz(-1.5900702) q[1];
rz(-0.20229907) q[3];
sx q[3];
rz(-1.9194366) q[3];
sx q[3];
rz(-1.5983456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8396987) q[2];
sx q[2];
rz(-1.7599186) q[2];
sx q[2];
rz(-1.1559486) q[2];
rz(0.65129748) q[3];
sx q[3];
rz(-2.7498701) q[3];
sx q[3];
rz(2.7294066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10525178) q[0];
sx q[0];
rz(-1.8008494) q[0];
sx q[0];
rz(-0.80129188) q[0];
rz(0.85805145) q[1];
sx q[1];
rz(-0.67887226) q[1];
sx q[1];
rz(-0.37989315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7004922) q[0];
sx q[0];
rz(-2.0394562) q[0];
sx q[0];
rz(0.20232135) q[0];
rz(-pi) q[1];
rz(0.11961898) q[2];
sx q[2];
rz(-1.3966171) q[2];
sx q[2];
rz(1.2816789) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26637883) q[1];
sx q[1];
rz(-0.76817487) q[1];
sx q[1];
rz(1.1298864) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44118446) q[3];
sx q[3];
rz(-1.5616186) q[3];
sx q[3];
rz(-1.6593195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6259367) q[2];
sx q[2];
rz(-2.8069324) q[2];
sx q[2];
rz(-0.77896172) q[2];
rz(2.7624779) q[3];
sx q[3];
rz(-1.4232676) q[3];
sx q[3];
rz(0.085722119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871386) q[0];
sx q[0];
rz(-1.4105281) q[0];
sx q[0];
rz(-2.6692303) q[0];
rz(-0.27913276) q[1];
sx q[1];
rz(-2.0886853) q[1];
sx q[1];
rz(-0.36769029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5147509) q[0];
sx q[0];
rz(-1.2184596) q[0];
sx q[0];
rz(-2.8295838) q[0];
rz(0.49007551) q[2];
sx q[2];
rz(-2.6948409) q[2];
sx q[2];
rz(1.2778145) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.20357) q[1];
sx q[1];
rz(-1.9782776) q[1];
sx q[1];
rz(-1.3111056) q[1];
rz(0.69766694) q[3];
sx q[3];
rz(-0.945745) q[3];
sx q[3];
rz(2.5755957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8184066) q[2];
sx q[2];
rz(-1.0419934) q[2];
sx q[2];
rz(-0.47427487) q[2];
rz(-0.91931528) q[3];
sx q[3];
rz(-1.5654469) q[3];
sx q[3];
rz(1.3924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11488386) q[0];
sx q[0];
rz(-1.3020496) q[0];
sx q[0];
rz(1.5337457) q[0];
rz(-2.907091) q[1];
sx q[1];
rz(-1.0954183) q[1];
sx q[1];
rz(-0.29874994) q[1];
rz(2.1420494) q[2];
sx q[2];
rz(-1.4951946) q[2];
sx q[2];
rz(-1.9517938) q[2];
rz(-0.017757105) q[3];
sx q[3];
rz(-1.5212301) q[3];
sx q[3];
rz(-2.7230079) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
