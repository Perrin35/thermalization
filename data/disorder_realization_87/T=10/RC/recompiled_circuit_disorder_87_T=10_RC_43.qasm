OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(5.5570931) q[0];
sx q[0];
rz(9.2232016) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(-0.5402686) q[1];
sx q[1];
rz(-0.93710605) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9286081) q[0];
sx q[0];
rz(-1.0092508) q[0];
sx q[0];
rz(-2.0748078) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0785157) q[2];
sx q[2];
rz(-2.1247851) q[2];
sx q[2];
rz(-0.89821494) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.571849) q[1];
sx q[1];
rz(-1.7535926) q[1];
sx q[1];
rz(0.079449541) q[1];
x q[2];
rz(1.3017544) q[3];
sx q[3];
rz(-1.6645414) q[3];
sx q[3];
rz(2.4274488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9822838) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(0.086159555) q[2];
rz(-0.75749767) q[3];
sx q[3];
rz(-2.3870654) q[3];
sx q[3];
rz(-2.0479726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(-2.5646097) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(-1.8151059) q[0];
rz(-1.2558698) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(-2.870141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46389929) q[0];
sx q[0];
rz(-1.9672183) q[0];
sx q[0];
rz(1.0248313) q[0];
rz(-pi) q[1];
rz(-1.2618622) q[2];
sx q[2];
rz(-1.075282) q[2];
sx q[2];
rz(2.3625284) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93310624) q[1];
sx q[1];
rz(-2.3098364) q[1];
sx q[1];
rz(2.0887124) q[1];
rz(-pi) q[2];
rz(-0.35630393) q[3];
sx q[3];
rz(-0.84173991) q[3];
sx q[3];
rz(-2.9475398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0469971) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(2.5644152) q[2];
rz(-2.2180637) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24401027) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(-2.999021) q[0];
rz(-1.7890731) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(-2.9325063) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81461834) q[0];
sx q[0];
rz(-2.1332392) q[0];
sx q[0];
rz(0.66280611) q[0];
rz(-pi) q[1];
rz(0.73968898) q[2];
sx q[2];
rz(-1.2581173) q[2];
sx q[2];
rz(-2.3198421) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0002034) q[1];
sx q[1];
rz(-2.5089426) q[1];
sx q[1];
rz(0.81694095) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61061065) q[3];
sx q[3];
rz(-2.7770677) q[3];
sx q[3];
rz(3.1162457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7769988) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(2.7374632) q[2];
rz(-1.2858307) q[3];
sx q[3];
rz(-1.1288246) q[3];
sx q[3];
rz(2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053112) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(0.96281111) q[0];
rz(2.6722233) q[1];
sx q[1];
rz(-2.5517187) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92361802) q[0];
sx q[0];
rz(-2.1305363) q[0];
sx q[0];
rz(-0.864242) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0175866) q[2];
sx q[2];
rz(-1.5713072) q[2];
sx q[2];
rz(-0.83204568) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5818134) q[1];
sx q[1];
rz(-2.1372791) q[1];
sx q[1];
rz(-0.15484667) q[1];
rz(-0.56941454) q[3];
sx q[3];
rz(-1.4953519) q[3];
sx q[3];
rz(0.30573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0756388) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(-0.11165079) q[2];
rz(0.81104898) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(-0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.951293) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(0.35650373) q[0];
rz(-2.6351392) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(-0.26062632) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8286752) q[0];
sx q[0];
rz(-0.14794359) q[0];
sx q[0];
rz(-0.28767985) q[0];
rz(-0.25116253) q[2];
sx q[2];
rz(-1.8166944) q[2];
sx q[2];
rz(-2.9159301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.645748) q[1];
sx q[1];
rz(-2.145088) q[1];
sx q[1];
rz(1.5917642) q[1];
x q[2];
rz(2.9959833) q[3];
sx q[3];
rz(-1.429261) q[3];
sx q[3];
rz(1.6573997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9115209) q[2];
sx q[2];
rz(-1.4804966) q[2];
sx q[2];
rz(-0.22932209) q[2];
rz(2.5991332) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(-2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77263537) q[0];
sx q[0];
rz(-2.2326523) q[0];
sx q[0];
rz(1.649958) q[0];
rz(-2.1024599) q[1];
sx q[1];
rz(-1.8455448) q[1];
sx q[1];
rz(1.414149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6468069) q[0];
sx q[0];
rz(-1.2959058) q[0];
sx q[0];
rz(-1.5216212) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8374149) q[2];
sx q[2];
rz(-1.068371) q[2];
sx q[2];
rz(-0.38602877) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.27281877) q[1];
sx q[1];
rz(-1.0827853) q[1];
sx q[1];
rz(-2.6150319) q[1];
rz(-3.066091) q[3];
sx q[3];
rz(-1.7864979) q[3];
sx q[3];
rz(-0.84738934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(-2.6489143) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(-1.3940575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3843) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(3.0199155) q[0];
rz(-1.1514459) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(-0.40245232) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1131176) q[0];
sx q[0];
rz(-2.3987781) q[0];
sx q[0];
rz(1.9755367) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0955663) q[2];
sx q[2];
rz(-1.5791025) q[2];
sx q[2];
rz(1.1549032) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9060045) q[1];
sx q[1];
rz(-1.8530122) q[1];
sx q[1];
rz(1.7547592) q[1];
rz(-pi) q[2];
rz(1.9190556) q[3];
sx q[3];
rz(-1.5446483) q[3];
sx q[3];
rz(2.1611283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.014331269) q[2];
sx q[2];
rz(-1.971259) q[2];
sx q[2];
rz(2.2793615) q[2];
rz(0.47752738) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(-0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35571337) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(-2.4086337) q[0];
rz(-0.14006242) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(-2.1070811) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5401934) q[0];
sx q[0];
rz(-1.7978151) q[0];
sx q[0];
rz(-2.9674203) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10266281) q[2];
sx q[2];
rz(-1.6409988) q[2];
sx q[2];
rz(-1.3753124) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1840399) q[1];
sx q[1];
rz(-1.5337481) q[1];
sx q[1];
rz(-2.1368105) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9018974) q[3];
sx q[3];
rz(-2.7125159) q[3];
sx q[3];
rz(-0.33340463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2404279) q[2];
sx q[2];
rz(-1.9463836) q[2];
sx q[2];
rz(-2.365716) q[2];
rz(0.72426978) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6999321) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(-3.124776) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.3341981) q[1];
sx q[1];
rz(-0.7787849) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6803857) q[0];
sx q[0];
rz(-1.9068423) q[0];
sx q[0];
rz(1.8340322) q[0];
rz(-1.0672827) q[2];
sx q[2];
rz(-2.6138966) q[2];
sx q[2];
rz(0.50349456) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15021819) q[1];
sx q[1];
rz(-2.7017936) q[1];
sx q[1];
rz(0.99779731) q[1];
rz(-pi) q[2];
rz(-3.0770244) q[3];
sx q[3];
rz(-1.350292) q[3];
sx q[3];
rz(-1.2793503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-2.4251535) q[2];
rz(1.6843494) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-0.15144908) q[0];
rz(0.17627136) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(0.72296468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6930981) q[0];
sx q[0];
rz(-1.4204331) q[0];
sx q[0];
rz(-0.10755121) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9777708) q[2];
sx q[2];
rz(-2.7343035) q[2];
sx q[2];
rz(-0.92330698) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55262676) q[1];
sx q[1];
rz(-0.50682658) q[1];
sx q[1];
rz(-1.6848906) q[1];
x q[2];
rz(2.74182) q[3];
sx q[3];
rz(-0.82953605) q[3];
sx q[3];
rz(-1.6403891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4621949) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(2.6160713) q[2];
rz(2.8578791) q[3];
sx q[3];
rz(-2.0339537) q[3];
sx q[3];
rz(-2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491966) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(-2.0422968) q[1];
sx q[1];
rz(-0.14930832) q[1];
sx q[1];
rz(2.2717448) q[1];
rz(0.47808403) q[2];
sx q[2];
rz(-2.1838084) q[2];
sx q[2];
rz(2.9391391) q[2];
rz(1.7952193) q[3];
sx q[3];
rz(-1.2524458) q[3];
sx q[3];
rz(2.8129775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
