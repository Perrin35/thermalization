OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(4.0868563) q[0];
sx q[0];
rz(9.950369) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(2.2367509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63431595) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(1.0213724) q[0];
rz(2.9154539) q[2];
sx q[2];
rz(-1.7625426) q[2];
sx q[2];
rz(-0.5118256) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2658087) q[1];
sx q[1];
rz(-1.8059397) q[1];
sx q[1];
rz(-1.9782515) q[1];
rz(-pi) q[2];
rz(1.8815133) q[3];
sx q[3];
rz(-1.7508535) q[3];
sx q[3];
rz(0.77375274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.20377542) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(3.0453483) q[2];
rz(-1.0359267) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7006943) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(-0.76517117) q[0];
rz(1.8493429) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(0.66295019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11766079) q[0];
sx q[0];
rz(-1.6370602) q[0];
sx q[0];
rz(-0.06225417) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0775954) q[2];
sx q[2];
rz(-0.71841824) q[2];
sx q[2];
rz(1.6990627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1124277) q[1];
sx q[1];
rz(-1.6532835) q[1];
sx q[1];
rz(0.44177456) q[1];
rz(-pi) q[2];
rz(1.7327659) q[3];
sx q[3];
rz(-2.0521945) q[3];
sx q[3];
rz(3.0961406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.092676) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(-0.66550231) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(0.22234017) q[0];
rz(-2.1242583) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(2.6229048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9496574) q[0];
sx q[0];
rz(-1.2149095) q[0];
sx q[0];
rz(-0.8215254) q[0];
rz(-pi) q[1];
rz(2.8027595) q[2];
sx q[2];
rz(-2.860184) q[2];
sx q[2];
rz(2.0237405) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7171214) q[1];
sx q[1];
rz(-2.3470504) q[1];
sx q[1];
rz(0.24309991) q[1];
rz(-pi) q[2];
x q[2];
rz(0.054611562) q[3];
sx q[3];
rz(-1.5715989) q[3];
sx q[3];
rz(-0.69566788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8609994) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-3.0730491) q[2];
rz(0.60244256) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(-0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.65072) q[0];
sx q[0];
rz(-0.77712178) q[0];
sx q[0];
rz(2.9673476) q[0];
rz(0.53025591) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(-0.51309103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78631567) q[0];
sx q[0];
rz(-2.516541) q[0];
sx q[0];
rz(3.0062208) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5330403) q[2];
sx q[2];
rz(-0.90494472) q[2];
sx q[2];
rz(-2.7666639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.40401134) q[1];
sx q[1];
rz(-2.0239132) q[1];
sx q[1];
rz(1.416942) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1107535) q[3];
sx q[3];
rz(-1.3645002) q[3];
sx q[3];
rz(-1.7145969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47485581) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(-0.22988698) q[2];
rz(0.41904467) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(-0.45927799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6722365) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(-2.4647734) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(-2.5255323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13000935) q[0];
sx q[0];
rz(-1.6066215) q[0];
sx q[0];
rz(-1.5183166) q[0];
x q[1];
rz(-1.6904171) q[2];
sx q[2];
rz(-1.9336485) q[2];
sx q[2];
rz(0.56064831) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7459813) q[1];
sx q[1];
rz(-0.99318722) q[1];
sx q[1];
rz(1.5773768) q[1];
rz(-pi) q[2];
rz(-1.945799) q[3];
sx q[3];
rz(-0.71205322) q[3];
sx q[3];
rz(2.860481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(2.8862254) q[2];
rz(1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(-0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(2.6690924) q[0];
rz(-2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(0.88476673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17381829) q[0];
sx q[0];
rz(-0.38823715) q[0];
sx q[0];
rz(1.3740747) q[0];
rz(-pi) q[1];
rz(-0.066263513) q[2];
sx q[2];
rz(-2.931086) q[2];
sx q[2];
rz(-1.9312242) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92237597) q[1];
sx q[1];
rz(-1.9799385) q[1];
sx q[1];
rz(-2.8273696) q[1];
rz(-pi) q[2];
rz(-1.0370449) q[3];
sx q[3];
rz(-1.3579206) q[3];
sx q[3];
rz(-0.45149976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3970967) q[2];
sx q[2];
rz(-0.8049736) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(1.7729676) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(3.080522) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(2.4087002) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7948579) q[0];
sx q[0];
rz(-1.1746527) q[0];
sx q[0];
rz(-2.5651155) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86874666) q[2];
sx q[2];
rz(-1.1940496) q[2];
sx q[2];
rz(1.4993315) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7072308) q[1];
sx q[1];
rz(-1.5773298) q[1];
sx q[1];
rz(2.1965501) q[1];
rz(-pi) q[2];
rz(-1.6781428) q[3];
sx q[3];
rz(-1.5296474) q[3];
sx q[3];
rz(2.6220208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0682893) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(2.627009) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.085389) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(-0.7094267) q[0];
rz(1.5052694) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(2.8628796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3364209) q[0];
sx q[0];
rz(-1.3609017) q[0];
sx q[0];
rz(-1.7091442) q[0];
x q[1];
rz(1.8144572) q[2];
sx q[2];
rz(-2.116034) q[2];
sx q[2];
rz(-1.3512163) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0875138) q[1];
sx q[1];
rz(-1.8808108) q[1];
sx q[1];
rz(2.9127496) q[1];
rz(-pi) q[2];
rz(1.7940984) q[3];
sx q[3];
rz(-1.6468862) q[3];
sx q[3];
rz(3.1267816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30148208) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(-3.0855132) q[2];
rz(0.85514832) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99496019) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(-2.1726998) q[0];
rz(0.47337198) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(-1.999058) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24507095) q[0];
sx q[0];
rz(-2.8327496) q[0];
sx q[0];
rz(1.8737428) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1387024) q[2];
sx q[2];
rz(-0.70334607) q[2];
sx q[2];
rz(-3.0386915) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.090230378) q[1];
sx q[1];
rz(-1.6962595) q[1];
sx q[1];
rz(-2.5999703) q[1];
rz(-pi) q[2];
rz(2.2409866) q[3];
sx q[3];
rz(-1.4108676) q[3];
sx q[3];
rz(-2.1538494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.896495) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(-1.0207821) q[2];
rz(2.8178689) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1616515) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(2.4401869) q[0];
rz(2.2258863) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(1.9030301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21047132) q[0];
sx q[0];
rz(-2.5744372) q[0];
sx q[0];
rz(-1.9803067) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6221223) q[2];
sx q[2];
rz(-0.81846279) q[2];
sx q[2];
rz(-2.3527956) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8348332) q[1];
sx q[1];
rz(-0.66654897) q[1];
sx q[1];
rz(-1.761761) q[1];
x q[2];
rz(1.402461) q[3];
sx q[3];
rz(-0.59909648) q[3];
sx q[3];
rz(-1.400791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68676585) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-0.043838538) q[2];
rz(1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(-2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.4713521) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(-0.025370601) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(-1.9699026) q[2];
sx q[2];
rz(-1.8532955) q[2];
sx q[2];
rz(2.4886139) q[2];
rz(0.79904859) q[3];
sx q[3];
rz(-1.6934762) q[3];
sx q[3];
rz(1.3140524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
