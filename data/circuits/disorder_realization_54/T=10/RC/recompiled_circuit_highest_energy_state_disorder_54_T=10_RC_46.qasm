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
rz(-0.71159166) q[0];
sx q[0];
rz(-2.0366259) q[0];
sx q[0];
rz(1.1397064) q[0];
rz(2.3501514) q[1];
sx q[1];
rz(-2.1005519) q[1];
sx q[1];
rz(-2.0120373) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68398841) q[0];
sx q[0];
rz(-2.2778371) q[0];
sx q[0];
rz(0.23914214) q[0];
rz(-pi) q[1];
rz(-1.5426251) q[2];
sx q[2];
rz(-0.90785387) q[2];
sx q[2];
rz(-2.4848795) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3870377) q[1];
sx q[1];
rz(-2.3554152) q[1];
sx q[1];
rz(-1.1494067) q[1];
rz(-pi) q[2];
rz(-0.33282354) q[3];
sx q[3];
rz(-1.6037919) q[3];
sx q[3];
rz(-0.28256114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1019885) q[2];
sx q[2];
rz(-0.81059376) q[2];
sx q[2];
rz(-2.2339036) q[2];
rz(3.0229819) q[3];
sx q[3];
rz(-1.6018931) q[3];
sx q[3];
rz(-2.8906726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949961) q[0];
sx q[0];
rz(-0.66221607) q[0];
sx q[0];
rz(-3.0372341) q[0];
rz(1.1907499) q[1];
sx q[1];
rz(-1.933681) q[1];
sx q[1];
rz(-2.6844535) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8123618) q[0];
sx q[0];
rz(-2.9456249) q[0];
sx q[0];
rz(-1.7295444) q[0];
rz(-2.1188583) q[2];
sx q[2];
rz(-0.75621683) q[2];
sx q[2];
rz(-0.05847419) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3981084) q[1];
sx q[1];
rz(-2.8116075) q[1];
sx q[1];
rz(-1.2742001) q[1];
rz(-pi) q[2];
rz(-2.2394286) q[3];
sx q[3];
rz(-1.4621825) q[3];
sx q[3];
rz(-1.4387552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5856058) q[2];
sx q[2];
rz(-1.1588691) q[2];
sx q[2];
rz(0.15698329) q[2];
rz(0.5213151) q[3];
sx q[3];
rz(-0.83955228) q[3];
sx q[3];
rz(0.014001525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3954725) q[0];
sx q[0];
rz(-0.59500256) q[0];
sx q[0];
rz(-0.34573063) q[0];
rz(-1.6861457) q[1];
sx q[1];
rz(-1.8691749) q[1];
sx q[1];
rz(1.570943) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60590505) q[0];
sx q[0];
rz(-0.55258026) q[0];
sx q[0];
rz(-1.7406169) q[0];
rz(-pi) q[1];
rz(0.94812553) q[2];
sx q[2];
rz(-1.5737163) q[2];
sx q[2];
rz(-0.44026431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6863043) q[1];
sx q[1];
rz(-1.6197816) q[1];
sx q[1];
rz(-2.9279885) q[1];
rz(-0.9216347) q[3];
sx q[3];
rz(-1.7329669) q[3];
sx q[3];
rz(2.0258486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3674783) q[2];
sx q[2];
rz(-2.9134637) q[2];
sx q[2];
rz(0.95743123) q[2];
rz(1.1227603) q[3];
sx q[3];
rz(-1.9072396) q[3];
sx q[3];
rz(1.3721589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289309) q[0];
sx q[0];
rz(-1.4224195) q[0];
sx q[0];
rz(-1.5976394) q[0];
rz(-1.3900025) q[1];
sx q[1];
rz(-1.3163687) q[1];
sx q[1];
rz(1.4768627) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2643971) q[0];
sx q[0];
rz(-2.0937284) q[0];
sx q[0];
rz(-1.8967129) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8091168) q[2];
sx q[2];
rz(-1.4891948) q[2];
sx q[2];
rz(-1.2343182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1424496) q[1];
sx q[1];
rz(-1.8987185) q[1];
sx q[1];
rz(-1.0161922) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8119726) q[3];
sx q[3];
rz(-1.8330036) q[3];
sx q[3];
rz(2.5548558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59820286) q[2];
sx q[2];
rz(-1.8566088) q[2];
sx q[2];
rz(-2.1447935) q[2];
rz(1.2654842) q[3];
sx q[3];
rz(-0.71804738) q[3];
sx q[3];
rz(-2.8640981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0604414) q[0];
sx q[0];
rz(-1.569898) q[0];
sx q[0];
rz(2.0695709) q[0];
rz(0.95442665) q[1];
sx q[1];
rz(-1.9642893) q[1];
sx q[1];
rz(1.6486453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.36906) q[0];
sx q[0];
rz(-2.323702) q[0];
sx q[0];
rz(-1.3007123) q[0];
rz(-pi) q[1];
rz(2.6482539) q[2];
sx q[2];
rz(-0.78560053) q[2];
sx q[2];
rz(2.0480905) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6449738) q[1];
sx q[1];
rz(-2.7826834) q[1];
sx q[1];
rz(2.7254057) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8478283) q[3];
sx q[3];
rz(-2.8076594) q[3];
sx q[3];
rz(1.6373843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2299049) q[2];
sx q[2];
rz(-2.7547084) q[2];
sx q[2];
rz(-1.9913199) q[2];
rz(-0.041042717) q[3];
sx q[3];
rz(-1.5951472) q[3];
sx q[3];
rz(-2.7400147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3732442) q[0];
sx q[0];
rz(-1.1078438) q[0];
sx q[0];
rz(1.1596229) q[0];
rz(-1.6805964) q[1];
sx q[1];
rz(-0.20080876) q[1];
sx q[1];
rz(-1.9083091) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.999812) q[0];
sx q[0];
rz(-1.2617869) q[0];
sx q[0];
rz(-2.3731493) q[0];
rz(-pi) q[1];
rz(0.75040011) q[2];
sx q[2];
rz(-3.0633492) q[2];
sx q[2];
rz(-0.83759826) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4602051) q[1];
sx q[1];
rz(-0.80599313) q[1];
sx q[1];
rz(-0.56350033) q[1];
x q[2];
rz(-1.308611) q[3];
sx q[3];
rz(-2.0823458) q[3];
sx q[3];
rz(1.2894323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.636574) q[2];
sx q[2];
rz(-2.7754112) q[2];
sx q[2];
rz(1.150307) q[2];
rz(-0.73927528) q[3];
sx q[3];
rz(-1.8940247) q[3];
sx q[3];
rz(0.44481835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66733852) q[0];
sx q[0];
rz(-0.69304729) q[0];
sx q[0];
rz(-0.44878238) q[0];
rz(-1.5397286) q[1];
sx q[1];
rz(-1.1069143) q[1];
sx q[1];
rz(-1.9805699) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13793531) q[0];
sx q[0];
rz(-1.9788392) q[0];
sx q[0];
rz(1.28027) q[0];
rz(-pi) q[1];
rz(2.6426605) q[2];
sx q[2];
rz(-2.009444) q[2];
sx q[2];
rz(-0.89887757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9185982) q[1];
sx q[1];
rz(-1.2410933) q[1];
sx q[1];
rz(-1.2597643) q[1];
rz(-0.26412873) q[3];
sx q[3];
rz(-1.0326516) q[3];
sx q[3];
rz(1.9028185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9728969) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(-1.9608344) q[2];
rz(2.0598038) q[3];
sx q[3];
rz(-1.6094145) q[3];
sx q[3];
rz(0.42657524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5931554) q[0];
sx q[0];
rz(-1.3189545) q[0];
sx q[0];
rz(0.68761188) q[0];
rz(1.9718735) q[1];
sx q[1];
rz(-2.1673188) q[1];
sx q[1];
rz(-2.7489472) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3340281) q[0];
sx q[0];
rz(-2.0453718) q[0];
sx q[0];
rz(-2.9827098) q[0];
rz(-pi) q[1];
rz(-2.7313825) q[2];
sx q[2];
rz(-0.33818076) q[2];
sx q[2];
rz(1.8800125) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6456522) q[1];
sx q[1];
rz(-2.0447746) q[1];
sx q[1];
rz(-1.4723634) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6527376) q[3];
sx q[3];
rz(-1.7590766) q[3];
sx q[3];
rz(0.57037607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16200599) q[2];
sx q[2];
rz(-2.1350828) q[2];
sx q[2];
rz(-1.4957734) q[2];
rz(1.6188072) q[3];
sx q[3];
rz(-2.3706172) q[3];
sx q[3];
rz(-0.60105598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9252121) q[0];
sx q[0];
rz(-2.7583211) q[0];
sx q[0];
rz(1.8076757) q[0];
rz(1.1434309) q[1];
sx q[1];
rz(-0.83088487) q[1];
sx q[1];
rz(-2.3480031) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.586602) q[0];
sx q[0];
rz(-1.7287917) q[0];
sx q[0];
rz(-1.4992856) q[0];
rz(-pi) q[1];
rz(1.4758395) q[2];
sx q[2];
rz(-0.61197241) q[2];
sx q[2];
rz(0.31699917) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55793437) q[1];
sx q[1];
rz(-2.0730632) q[1];
sx q[1];
rz(1.7947547) q[1];
rz(-1.3511485) q[3];
sx q[3];
rz(-2.6815412) q[3];
sx q[3];
rz(-1.8203514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6952343) q[2];
sx q[2];
rz(-2.0909205) q[2];
sx q[2];
rz(1.0515155) q[2];
rz(-2.4255883) q[3];
sx q[3];
rz(-0.99596888) q[3];
sx q[3];
rz(-0.21014617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150772) q[0];
sx q[0];
rz(-2.2861013) q[0];
sx q[0];
rz(-0.18095428) q[0];
rz(1.1894233) q[1];
sx q[1];
rz(-1.1782497) q[1];
sx q[1];
rz(0.98035556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086209379) q[0];
sx q[0];
rz(-0.50859857) q[0];
sx q[0];
rz(-0.17091708) q[0];
rz(-pi) q[1];
rz(-0.75677121) q[2];
sx q[2];
rz(-0.59322651) q[2];
sx q[2];
rz(2.7293918) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7804523) q[1];
sx q[1];
rz(-1.3899511) q[1];
sx q[1];
rz(-2.3922582) q[1];
x q[2];
rz(0.47360955) q[3];
sx q[3];
rz(-2.7033812) q[3];
sx q[3];
rz(-1.1058863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1159726) q[2];
sx q[2];
rz(-0.16416922) q[2];
sx q[2];
rz(1.8170961) q[2];
rz(-0.65274158) q[3];
sx q[3];
rz(-0.96941152) q[3];
sx q[3];
rz(-3.1033707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7272335) q[0];
sx q[0];
rz(-1.7200732) q[0];
sx q[0];
rz(2.1678069) q[0];
rz(-2.4173792) q[1];
sx q[1];
rz(-1.0754633) q[1];
sx q[1];
rz(0.16574688) q[1];
rz(-1.235511) q[2];
sx q[2];
rz(-1.4852471) q[2];
sx q[2];
rz(0.5813364) q[2];
rz(-0.11304819) q[3];
sx q[3];
rz(-0.77533508) q[3];
sx q[3];
rz(-1.2067457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
