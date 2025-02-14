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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6467616) q[0];
sx q[0];
rz(-1.1755623) q[0];
sx q[0];
rz(0.075972164) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5025592) q[2];
sx q[2];
rz(-2.0026853) q[2];
sx q[2];
rz(-2.2186861) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9153629) q[1];
sx q[1];
rz(-2.1945476) q[1];
sx q[1];
rz(2.8585494) q[1];
rz(-pi) q[2];
rz(1.3320311) q[3];
sx q[3];
rz(-0.70203915) q[3];
sx q[3];
rz(0.89621893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4932058) q[2];
sx q[2];
rz(-1.9923261) q[2];
sx q[2];
rz(-3.1374078) q[2];
rz(-0.74797136) q[3];
sx q[3];
rz(-2.3145521) q[3];
sx q[3];
rz(-2.9048257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0325322) q[0];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4499265) q[0];
sx q[0];
rz(-2.0332893) q[0];
sx q[0];
rz(1.7989547) q[0];
x q[1];
rz(0.02702464) q[2];
sx q[2];
rz(-1.4725497) q[2];
sx q[2];
rz(-0.25729986) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8557094) q[1];
sx q[1];
rz(-2.6220589) q[1];
sx q[1];
rz(2.829771) q[1];
rz(-pi) q[2];
rz(-0.22547884) q[3];
sx q[3];
rz(-1.6957458) q[3];
sx q[3];
rz(-1.8285816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8740497) q[2];
sx q[2];
rz(-2.0127313) q[2];
sx q[2];
rz(1.4675325) q[2];
rz(-2.2360146) q[3];
sx q[3];
rz(-1.0970683) q[3];
sx q[3];
rz(2.9300698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8242789) q[0];
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
rz(1.9267265) q[0];
sx q[0];
rz(-1.1385554) q[0];
sx q[0];
rz(1.3316447) q[0];
rz(0.27806313) q[2];
sx q[2];
rz(-0.5812656) q[2];
sx q[2];
rz(-0.9424302) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7819848) q[1];
sx q[1];
rz(-2.4949412) q[1];
sx q[1];
rz(-2/(3*pi)) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0140225) q[3];
sx q[3];
rz(-2.2446998) q[3];
sx q[3];
rz(1.2942838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5182284) q[2];
sx q[2];
rz(-0.48064226) q[2];
sx q[2];
rz(2.9474958) q[2];
rz(2.569681) q[3];
sx q[3];
rz(-1.5827936) q[3];
sx q[3];
rz(-1.5552049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67683515) q[0];
sx q[0];
rz(-2.0763626) q[0];
sx q[0];
rz(-1.42365) q[0];
rz(2.8061197) q[1];
sx q[1];
rz(-0.60233855) q[1];
sx q[1];
rz(0.64340341) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5455824) q[0];
sx q[0];
rz(-1.032858) q[0];
sx q[0];
rz(0.62278231) q[0];
rz(0.14428987) q[2];
sx q[2];
rz(-0.57960287) q[2];
sx q[2];
rz(1.6336827) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57694084) q[1];
sx q[1];
rz(-1.1168861) q[1];
sx q[1];
rz(0.79344016) q[1];
rz(-pi) q[2];
rz(2.0213177) q[3];
sx q[3];
rz(-1.4304597) q[3];
sx q[3];
rz(-0.72718638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.246835) q[2];
sx q[2];
rz(-1.094123) q[2];
sx q[2];
rz(-0.82478729) q[2];
rz(-2.5347533) q[3];
sx q[3];
rz(-1.2068799) q[3];
sx q[3];
rz(1.3185893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7428699) q[0];
sx q[0];
rz(-0.56951183) q[0];
sx q[0];
rz(-0.2670162) q[0];
rz(0.46785242) q[1];
sx q[1];
rz(-2.0303969) q[1];
sx q[1];
rz(1.1763447) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6716135) q[0];
sx q[0];
rz(-1.592127) q[0];
sx q[0];
rz(1.6213378) q[0];
rz(2.3651334) q[2];
sx q[2];
rz(-1.0506223) q[2];
sx q[2];
rz(0.3683683) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49217192) q[1];
sx q[1];
rz(-0.4539868) q[1];
sx q[1];
rz(1.4439871) q[1];
rz(0.98988804) q[3];
sx q[3];
rz(-2.2047289) q[3];
sx q[3];
rz(0.1617728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0331369) q[2];
sx q[2];
rz(-1.5908073) q[2];
sx q[2];
rz(-0.94672686) q[2];
rz(1.4416134) q[3];
sx q[3];
rz(-1.0333034) q[3];
sx q[3];
rz(-1.3942963) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.5627728) q[1];
sx q[1];
rz(1.2616166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3917106) q[0];
sx q[0];
rz(-0.80750033) q[0];
sx q[0];
rz(0.9901643) q[0];
rz(-pi) q[1];
rz(-0.67503898) q[2];
sx q[2];
rz(-1.0770105) q[2];
sx q[2];
rz(-1.5641649) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4289866) q[1];
sx q[1];
rz(-1.3416051) q[1];
sx q[1];
rz(1.5043753) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2618939) q[3];
sx q[3];
rz(-2.8692103) q[3];
sx q[3];
rz(0.72379823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6235846) q[2];
sx q[2];
rz(-1.1058747) q[2];
sx q[2];
rz(2.8688431) q[2];
rz(0.92877156) q[3];
sx q[3];
rz(-2.2835968) q[3];
sx q[3];
rz(-1.4906918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6269094) q[0];
sx q[0];
rz(-1.0310443) q[0];
sx q[0];
rz(-2.2169901) q[0];
rz(-0.54479105) q[1];
sx q[1];
rz(-1.7926615) q[1];
sx q[1];
rz(-0.10442385) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4686779) q[0];
sx q[0];
rz(-0.12355655) q[0];
sx q[0];
rz(0.79609032) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41506501) q[2];
sx q[2];
rz(-1.7156148) q[2];
sx q[2];
rz(0.96127779) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.663781) q[1];
sx q[1];
rz(-0.38214499) q[1];
sx q[1];
rz(-0.93375979) q[1];
rz(-pi) q[2];
rz(0.35108836) q[3];
sx q[3];
rz(-1.2005873) q[3];
sx q[3];
rz(0.58634392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32120785) q[2];
sx q[2];
rz(-1.2106004) q[2];
sx q[2];
rz(-1.661181) q[2];
rz(0.004247578) q[3];
sx q[3];
rz(-2.28528) q[3];
sx q[3];
rz(-0.64928865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55411196) q[0];
sx q[0];
rz(-0.09621796) q[0];
sx q[0];
rz(-1.8884678) q[0];
rz(-2.9907277) q[1];
sx q[1];
rz(-0.8395218) q[1];
sx q[1];
rz(-1.2996659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2723994) q[0];
sx q[0];
rz(-1.0080833) q[0];
sx q[0];
rz(1.2666232) q[0];
x q[1];
rz(-1.2339044) q[2];
sx q[2];
rz(-0.46122069) q[2];
sx q[2];
rz(2.5835844) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7598301) q[1];
sx q[1];
rz(-0.88517939) q[1];
sx q[1];
rz(1.5515224) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9392936) q[3];
sx q[3];
rz(-1.9194366) q[3];
sx q[3];
rz(-1.5983456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3018939) q[2];
sx q[2];
rz(-1.7599186) q[2];
sx q[2];
rz(1.1559486) q[2];
rz(-2.4902952) q[3];
sx q[3];
rz(-0.39172253) q[3];
sx q[3];
rz(0.41218606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0363409) q[0];
sx q[0];
rz(-1.8008494) q[0];
sx q[0];
rz(0.80129188) q[0];
rz(-2.2835412) q[1];
sx q[1];
rz(-2.4627204) q[1];
sx q[1];
rz(-2.7616995) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86738619) q[0];
sx q[0];
rz(-0.50747061) q[0];
sx q[0];
rz(1.9486289) q[0];
rz(-1.7462037) q[2];
sx q[2];
rz(-1.6885969) q[2];
sx q[2];
rz(-2.8316488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26637883) q[1];
sx q[1];
rz(-2.3734178) q[1];
sx q[1];
rz(-1.1298864) q[1];
rz(-pi) q[2];
rz(-1.5809459) q[3];
sx q[3];
rz(-2.0119609) q[3];
sx q[3];
rz(-0.084189296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6259367) q[2];
sx q[2];
rz(-2.8069324) q[2];
sx q[2];
rz(-0.77896172) q[2];
rz(2.7624779) q[3];
sx q[3];
rz(-1.718325) q[3];
sx q[3];
rz(3.0558705) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15445408) q[0];
sx q[0];
rz(-1.7310646) q[0];
sx q[0];
rz(2.6692303) q[0];
rz(-0.27913276) q[1];
sx q[1];
rz(-2.0886853) q[1];
sx q[1];
rz(2.7739024) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.37848) q[0];
sx q[0];
rz(-2.6753278) q[0];
sx q[0];
rz(0.87514042) q[0];
rz(-pi) q[1];
rz(-1.7925749) q[2];
sx q[2];
rz(-1.9618755) q[2];
sx q[2];
rz(-1.329601) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34741857) q[1];
sx q[1];
rz(-2.6623335) q[1];
sx q[1];
rz(2.6049331) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69766694) q[3];
sx q[3];
rz(-0.945745) q[3];
sx q[3];
rz(0.56599697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3231861) q[2];
sx q[2];
rz(-1.0419934) q[2];
sx q[2];
rz(-0.47427487) q[2];
rz(0.91931528) q[3];
sx q[3];
rz(-1.5654469) q[3];
sx q[3];
rz(1.7491755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267088) q[0];
sx q[0];
rz(-1.8395431) q[0];
sx q[0];
rz(-1.607847) q[0];
rz(-2.907091) q[1];
sx q[1];
rz(-1.0954183) q[1];
sx q[1];
rz(-0.29874994) q[1];
rz(-2.1420494) q[2];
sx q[2];
rz(-1.6463981) q[2];
sx q[2];
rz(1.1897988) q[2];
rz(0.017757105) q[3];
sx q[3];
rz(-1.6203625) q[3];
sx q[3];
rz(0.41858471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
