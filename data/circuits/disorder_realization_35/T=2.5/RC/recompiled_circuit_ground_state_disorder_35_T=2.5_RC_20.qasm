OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93249455) q[0];
sx q[0];
rz(-0.068109186) q[0];
sx q[0];
rz(0.77686247) q[0];
rz(-1.3070973) q[1];
sx q[1];
rz(-0.96944648) q[1];
sx q[1];
rz(1.3219272) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.043524) q[0];
sx q[0];
rz(-0.72854818) q[0];
sx q[0];
rz(-1.3619955) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.916531) q[2];
sx q[2];
rz(-2.0421197) q[2];
sx q[2];
rz(2.2668348) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2474587) q[1];
sx q[1];
rz(-1.1620635) q[1];
sx q[1];
rz(3.1097059) q[1];
rz(0.98422756) q[3];
sx q[3];
rz(-2.7968458) q[3];
sx q[3];
rz(-1.00351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.063735828) q[2];
sx q[2];
rz(-1.3602164) q[2];
sx q[2];
rz(-2.784101) q[2];
rz(2.5743971) q[3];
sx q[3];
rz(-1.4128069) q[3];
sx q[3];
rz(0.43958694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95328632) q[0];
sx q[0];
rz(-2.8585241) q[0];
sx q[0];
rz(1.333492) q[0];
rz(0.4370583) q[1];
sx q[1];
rz(-0.60085618) q[1];
sx q[1];
rz(2.3016047) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8020426) q[0];
sx q[0];
rz(-1.4176482) q[0];
sx q[0];
rz(1.0569554) q[0];
rz(-pi) q[1];
rz(0.59462585) q[2];
sx q[2];
rz(-1.42808) q[2];
sx q[2];
rz(-3.0171301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6691554) q[1];
sx q[1];
rz(-1.0095114) q[1];
sx q[1];
rz(-0.88772933) q[1];
x q[2];
rz(1.6343104) q[3];
sx q[3];
rz(-0.5507142) q[3];
sx q[3];
rz(2.7515819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4240894) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(-1.4060414) q[2];
rz(2.9794335) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(2.6774008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9185987) q[0];
sx q[0];
rz(-3.0610237) q[0];
sx q[0];
rz(0.42174569) q[0];
rz(0.84284198) q[1];
sx q[1];
rz(-1.3858567) q[1];
sx q[1];
rz(-0.27580321) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0378135) q[0];
sx q[0];
rz(-2.7386754) q[0];
sx q[0];
rz(2.518464) q[0];
x q[1];
rz(-2.7245443) q[2];
sx q[2];
rz(-0.50784238) q[2];
sx q[2];
rz(-1.6091953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0581857) q[1];
sx q[1];
rz(-0.32489932) q[1];
sx q[1];
rz(-0.16375457) q[1];
x q[2];
rz(0.86902501) q[3];
sx q[3];
rz(-1.4956253) q[3];
sx q[3];
rz(1.6628671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3383823) q[2];
sx q[2];
rz(-2.3440177) q[2];
sx q[2];
rz(2.4513643) q[2];
rz(2.7817182) q[3];
sx q[3];
rz(-1.7706324) q[3];
sx q[3];
rz(-2.7580822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.2279219) q[0];
sx q[0];
rz(-2.0991195) q[0];
sx q[0];
rz(-0.87164718) q[0];
rz(2.7323885) q[1];
sx q[1];
rz(-0.49574655) q[1];
sx q[1];
rz(-0.31235487) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.384639) q[0];
sx q[0];
rz(-2.0142382) q[0];
sx q[0];
rz(-1.7009981) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6705475) q[2];
sx q[2];
rz(-1.776223) q[2];
sx q[2];
rz(-2.2769986) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96906584) q[1];
sx q[1];
rz(-1.2102264) q[1];
sx q[1];
rz(1.3532624) q[1];
rz(-pi) q[2];
rz(2.4931025) q[3];
sx q[3];
rz(-1.7269508) q[3];
sx q[3];
rz(-2.3462669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.006134) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(0.82497605) q[2];
rz(-1.4247591) q[3];
sx q[3];
rz(-2.8202839) q[3];
sx q[3];
rz(2.7062374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.41998) q[0];
sx q[0];
rz(-1.551832) q[0];
sx q[0];
rz(-2.2671674) q[0];
rz(-0.32886109) q[1];
sx q[1];
rz(-1.3855653) q[1];
sx q[1];
rz(-2.0430476) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4916353) q[0];
sx q[0];
rz(-1.1960317) q[0];
sx q[0];
rz(1.4521506) q[0];
rz(-pi) q[1];
rz(-0.036089049) q[2];
sx q[2];
rz(-0.83613013) q[2];
sx q[2];
rz(-3.0826838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9096127) q[1];
sx q[1];
rz(-1.8706338) q[1];
sx q[1];
rz(-1.058387) q[1];
rz(-pi) q[2];
rz(-0.9807113) q[3];
sx q[3];
rz(-1.8030589) q[3];
sx q[3];
rz(1.6764318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.141779) q[2];
sx q[2];
rz(-1.9580611) q[2];
sx q[2];
rz(-0.67406526) q[2];
rz(1.6541803) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(-0.28356799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0312626) q[0];
sx q[0];
rz(-2.1188348) q[0];
sx q[0];
rz(-1.7559825) q[0];
rz(-1.1194057) q[1];
sx q[1];
rz(-1.6740572) q[1];
sx q[1];
rz(-1.7562235) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35753174) q[0];
sx q[0];
rz(-0.94793944) q[0];
sx q[0];
rz(2.8297367) q[0];
rz(1.4848861) q[2];
sx q[2];
rz(-1.0738147) q[2];
sx q[2];
rz(-2.6045406) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5934903) q[1];
sx q[1];
rz(-1.8917276) q[1];
sx q[1];
rz(-1.7744508) q[1];
x q[2];
rz(-2.9117111) q[3];
sx q[3];
rz(-0.15860367) q[3];
sx q[3];
rz(-1.6386709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.063227) q[2];
sx q[2];
rz(-0.93457064) q[2];
sx q[2];
rz(-2.7654977) q[2];
rz(2.0070576) q[3];
sx q[3];
rz(-2.7781656) q[3];
sx q[3];
rz(0.80662066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.019808708) q[0];
sx q[0];
rz(-0.60096318) q[0];
sx q[0];
rz(-3.0390749) q[0];
rz(-2.5566697) q[1];
sx q[1];
rz(-2.1939317) q[1];
sx q[1];
rz(1.1835416) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7022423) q[0];
sx q[0];
rz(-2.9493414) q[0];
sx q[0];
rz(-0.73701145) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7223515) q[2];
sx q[2];
rz(-0.70690522) q[2];
sx q[2];
rz(-2.8173994) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.592748) q[1];
sx q[1];
rz(-1.8251437) q[1];
sx q[1];
rz(-2.2266822) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3041586) q[3];
sx q[3];
rz(-2.023306) q[3];
sx q[3];
rz(-3.0686015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3465603) q[2];
sx q[2];
rz(-2.2117895) q[2];
sx q[2];
rz(-2.9417876) q[2];
rz(-2.0810769) q[3];
sx q[3];
rz(-2.6001055) q[3];
sx q[3];
rz(-0.04960355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26315966) q[0];
sx q[0];
rz(-1.6596154) q[0];
sx q[0];
rz(-0.31295452) q[0];
rz(-0.9494268) q[1];
sx q[1];
rz(-1.3525617) q[1];
sx q[1];
rz(-1.4535646) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4623186) q[0];
sx q[0];
rz(-1.3841713) q[0];
sx q[0];
rz(-1.4786722) q[0];
rz(-1.4916999) q[2];
sx q[2];
rz(-2.1888615) q[2];
sx q[2];
rz(-1.5072418) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9447228) q[1];
sx q[1];
rz(-2.2579282) q[1];
sx q[1];
rz(1.2792759) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9096776) q[3];
sx q[3];
rz(-1.1383071) q[3];
sx q[3];
rz(-2.5755455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35385418) q[2];
sx q[2];
rz(-0.9757897) q[2];
sx q[2];
rz(-1.806951) q[2];
rz(-1.8169962) q[3];
sx q[3];
rz(-2.4748804) q[3];
sx q[3];
rz(-1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1052833) q[0];
sx q[0];
rz(-0.61573017) q[0];
sx q[0];
rz(0.69806725) q[0];
rz(2.2659194) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(2.7811513) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7822736) q[0];
sx q[0];
rz(-2.3397315) q[0];
sx q[0];
rz(0.62863056) q[0];
rz(-2.9827732) q[2];
sx q[2];
rz(-0.38714001) q[2];
sx q[2];
rz(0.57986812) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7061758) q[1];
sx q[1];
rz(-0.98505322) q[1];
sx q[1];
rz(2.6800568) q[1];
x q[2];
rz(-1.4029986) q[3];
sx q[3];
rz(-2.1133964) q[3];
sx q[3];
rz(-2.8436273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2299819) q[2];
sx q[2];
rz(-1.4260099) q[2];
sx q[2];
rz(-0.81725517) q[2];
rz(-2.3115555) q[3];
sx q[3];
rz(-0.54098141) q[3];
sx q[3];
rz(-0.219492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5068186) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(-3.1095374) q[0];
rz(1.8219148) q[1];
sx q[1];
rz(-1.7664884) q[1];
sx q[1];
rz(2.4874036) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35857707) q[0];
sx q[0];
rz(-2.9111283) q[0];
sx q[0];
rz(2.5925267) q[0];
x q[1];
rz(1.1172574) q[2];
sx q[2];
rz(-0.83570601) q[2];
sx q[2];
rz(2.105956) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.75958222) q[1];
sx q[1];
rz(-2.0339662) q[1];
sx q[1];
rz(2.5968371) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8579162) q[3];
sx q[3];
rz(-0.95796889) q[3];
sx q[3];
rz(-0.046663849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0111982) q[2];
sx q[2];
rz(-1.0426714) q[2];
sx q[2];
rz(-2.1350258) q[2];
rz(0.69342962) q[3];
sx q[3];
rz(-1.8208241) q[3];
sx q[3];
rz(-1.2815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4257767) q[0];
sx q[0];
rz(-0.67831138) q[0];
sx q[0];
rz(-0.074445733) q[0];
rz(0.88059942) q[1];
sx q[1];
rz(-2.0279299) q[1];
sx q[1];
rz(0.50150064) q[1];
rz(-2.1660317) q[2];
sx q[2];
rz(-0.79339334) q[2];
sx q[2];
rz(-1.7421772) q[2];
rz(2.4453658) q[3];
sx q[3];
rz(-1.8693557) q[3];
sx q[3];
rz(-2.1341677) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
