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
rz(-2.969279) q[0];
sx q[0];
rz(-0.20792374) q[0];
sx q[0];
rz(1.8508258) q[0];
rz(0.15200226) q[1];
sx q[1];
rz(3.7842964) q[1];
sx q[1];
rz(9.8531129) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9243619) q[0];
sx q[0];
rz(-0.15146449) q[0];
sx q[0];
rz(2.3343655) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4691888) q[2];
sx q[2];
rz(-1.8032296) q[2];
sx q[2];
rz(-1.642579) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2178035) q[1];
sx q[1];
rz(-2.3084967) q[1];
sx q[1];
rz(-1.6290725) q[1];
x q[2];
rz(1.3702014) q[3];
sx q[3];
rz(-1.7852968) q[3];
sx q[3];
rz(0.89591276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20994818) q[2];
sx q[2];
rz(-2.1817709) q[2];
sx q[2];
rz(1.9317365) q[2];
rz(-0.079553902) q[3];
sx q[3];
rz(-0.39593655) q[3];
sx q[3];
rz(2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2570268) q[0];
sx q[0];
rz(-0.50978065) q[0];
sx q[0];
rz(0.38401815) q[0];
rz(2.2513023) q[1];
sx q[1];
rz(-1.9046116) q[1];
sx q[1];
rz(-2.8382137) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10643364) q[0];
sx q[0];
rz(-2.1037397) q[0];
sx q[0];
rz(-0.22179166) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0069507) q[2];
sx q[2];
rz(-0.56098962) q[2];
sx q[2];
rz(-0.84535384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0413826) q[1];
sx q[1];
rz(-0.1100556) q[1];
sx q[1];
rz(-0.96918126) q[1];
rz(1.1214662) q[3];
sx q[3];
rz(-2.1746965) q[3];
sx q[3];
rz(0.7627129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7552135) q[2];
sx q[2];
rz(-1.1073802) q[2];
sx q[2];
rz(0.73491043) q[2];
rz(0.84878659) q[3];
sx q[3];
rz(-2.3990192) q[3];
sx q[3];
rz(-0.36062226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81247771) q[0];
sx q[0];
rz(-2.766093) q[0];
sx q[0];
rz(3.062881) q[0];
rz(-2.3178237) q[1];
sx q[1];
rz(-2.0562101) q[1];
sx q[1];
rz(-1.2944006) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1549653) q[0];
sx q[0];
rz(-1.6324348) q[0];
sx q[0];
rz(-3.0810647) q[0];
rz(-0.95125385) q[2];
sx q[2];
rz(-1.9980556) q[2];
sx q[2];
rz(-2.6841629) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30101206) q[1];
sx q[1];
rz(-2.8361179) q[1];
sx q[1];
rz(-0.39922797) q[1];
x q[2];
rz(-1.8494959) q[3];
sx q[3];
rz(-0.72360814) q[3];
sx q[3];
rz(2.0268745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2417629) q[2];
sx q[2];
rz(-0.37909847) q[2];
sx q[2];
rz(2.3251593) q[2];
rz(1.624931) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(2.4412156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65275943) q[0];
sx q[0];
rz(-0.82010287) q[0];
sx q[0];
rz(0.56831992) q[0];
rz(-1.9991416) q[1];
sx q[1];
rz(-1.3520974) q[1];
sx q[1];
rz(1.6894587) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84684337) q[0];
sx q[0];
rz(-2.1921625) q[0];
sx q[0];
rz(-1.2970379) q[0];
rz(-2.9617519) q[2];
sx q[2];
rz(-1.4786063) q[2];
sx q[2];
rz(-1.8123019) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.3209705) q[1];
sx q[1];
rz(-1.8838716) q[1];
sx q[1];
rz(0.93108196) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0978974) q[3];
sx q[3];
rz(-2.7644025) q[3];
sx q[3];
rz(0.24172129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51231724) q[2];
sx q[2];
rz(-2.6111111) q[2];
sx q[2];
rz(-3.0647035) q[2];
rz(-2.7422089) q[3];
sx q[3];
rz(-0.9136343) q[3];
sx q[3];
rz(-2.5975749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15453108) q[0];
sx q[0];
rz(-0.35683826) q[0];
sx q[0];
rz(-0.29990184) q[0];
rz(1.1497644) q[1];
sx q[1];
rz(-2.1297784) q[1];
sx q[1];
rz(0.48847517) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91575275) q[0];
sx q[0];
rz(-1.7560648) q[0];
sx q[0];
rz(3.1082694) q[0];
rz(-0.68814338) q[2];
sx q[2];
rz(-0.71758413) q[2];
sx q[2];
rz(-2.6731499) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9087968) q[1];
sx q[1];
rz(-0.98728647) q[1];
sx q[1];
rz(1.5481434) q[1];
rz(-2.4998922) q[3];
sx q[3];
rz(-1.3535168) q[3];
sx q[3];
rz(-0.02441306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38146314) q[2];
sx q[2];
rz(-3.0035786) q[2];
sx q[2];
rz(-1.6628954) q[2];
rz(2.0315157) q[3];
sx q[3];
rz(-0.9980945) q[3];
sx q[3];
rz(-2.4983675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4019302) q[0];
sx q[0];
rz(-1.1818385) q[0];
sx q[0];
rz(0.31841835) q[0];
rz(-0.034612522) q[1];
sx q[1];
rz(-2.5739659) q[1];
sx q[1];
rz(2.6908223) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50758963) q[0];
sx q[0];
rz(-1.5566155) q[0];
sx q[0];
rz(2.0812278) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87574739) q[2];
sx q[2];
rz(-1.0330457) q[2];
sx q[2];
rz(1.8132096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.602421) q[1];
sx q[1];
rz(-0.97214593) q[1];
sx q[1];
rz(-2.6461698) q[1];
x q[2];
rz(1.9616021) q[3];
sx q[3];
rz(-2.3473661) q[3];
sx q[3];
rz(-0.022217928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.047711756) q[2];
sx q[2];
rz(-0.1897976) q[2];
sx q[2];
rz(3.0799358) q[2];
rz(-0.098585248) q[3];
sx q[3];
rz(-0.74461377) q[3];
sx q[3];
rz(-0.70257598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3848569) q[0];
sx q[0];
rz(-2.0282133) q[0];
sx q[0];
rz(0.039948832) q[0];
rz(-0.7705676) q[1];
sx q[1];
rz(-2.4295085) q[1];
sx q[1];
rz(-0.22824731) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6421445) q[0];
sx q[0];
rz(-0.11760437) q[0];
sx q[0];
rz(2.2884877) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8769774) q[2];
sx q[2];
rz(-1.5723229) q[2];
sx q[2];
rz(2.7806139) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.764594) q[1];
sx q[1];
rz(-1.6239163) q[1];
sx q[1];
rz(-3.0594917) q[1];
rz(-pi) q[2];
rz(2.7258148) q[3];
sx q[3];
rz(-2.1630187) q[3];
sx q[3];
rz(-1.3982738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.063529) q[2];
sx q[2];
rz(-2.0827796) q[2];
sx q[2];
rz(0.25827363) q[2];
rz(-2.9410948) q[3];
sx q[3];
rz(-0.79810464) q[3];
sx q[3];
rz(0.18331461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9789326) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(-2.8343416) q[0];
rz(-1.2112674) q[1];
sx q[1];
rz(-0.4937506) q[1];
sx q[1];
rz(0.58120751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4330313) q[0];
sx q[0];
rz(-1.975346) q[0];
sx q[0];
rz(-0.033971196) q[0];
x q[1];
rz(1.4072335) q[2];
sx q[2];
rz(-0.86771655) q[2];
sx q[2];
rz(2.475098) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8381184) q[1];
sx q[1];
rz(-1.7707157) q[1];
sx q[1];
rz(0.78472991) q[1];
rz(-pi) q[2];
rz(2.5804306) q[3];
sx q[3];
rz(-2.5303839) q[3];
sx q[3];
rz(-1.2964013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6792949) q[2];
sx q[2];
rz(-2.0079948) q[2];
sx q[2];
rz(-0.031166859) q[2];
rz(2.9160685) q[3];
sx q[3];
rz(-1.3616819) q[3];
sx q[3];
rz(-1.0303191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0533326) q[0];
sx q[0];
rz(-0.044476155) q[0];
sx q[0];
rz(2.4326676) q[0];
rz(-2.9290579) q[1];
sx q[1];
rz(-1.0147076) q[1];
sx q[1];
rz(0.85740352) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5564726) q[0];
sx q[0];
rz(-1.5758638) q[0];
sx q[0];
rz(-1.7005672) q[0];
rz(2.6231758) q[2];
sx q[2];
rz(-0.5974434) q[2];
sx q[2];
rz(0.13472508) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9136476) q[1];
sx q[1];
rz(-1.040732) q[1];
sx q[1];
rz(1.1157406) q[1];
rz(-0.69714947) q[3];
sx q[3];
rz(-1.7350041) q[3];
sx q[3];
rz(2.177161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7015486) q[2];
sx q[2];
rz(-2.7893119) q[2];
sx q[2];
rz(0.80553833) q[2];
rz(-2.7696179) q[3];
sx q[3];
rz(-1.6476846) q[3];
sx q[3];
rz(0.39723799) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1086248) q[0];
sx q[0];
rz(-1.2297577) q[0];
sx q[0];
rz(-2.1627872) q[0];
rz(0.72273123) q[1];
sx q[1];
rz(-2.0131854) q[1];
sx q[1];
rz(0.61789787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.582983) q[0];
sx q[0];
rz(-1.8454058) q[0];
sx q[0];
rz(0.78297575) q[0];
rz(-pi) q[1];
rz(0.65872569) q[2];
sx q[2];
rz(-1.9945696) q[2];
sx q[2];
rz(1.2032594) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56983518) q[1];
sx q[1];
rz(-1.7392842) q[1];
sx q[1];
rz(-1.7525826) q[1];
x q[2];
rz(-0.79094751) q[3];
sx q[3];
rz(-0.11868782) q[3];
sx q[3];
rz(3.0968248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.939398) q[2];
sx q[2];
rz(-1.3556182) q[2];
sx q[2];
rz(0.00016577684) q[2];
rz(-0.58445066) q[3];
sx q[3];
rz(-1.0060468) q[3];
sx q[3];
rz(0.59529006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7547739) q[0];
sx q[0];
rz(-1.5703572) q[0];
sx q[0];
rz(1.5686709) q[0];
rz(1.3407002) q[1];
sx q[1];
rz(-2.0410213) q[1];
sx q[1];
rz(-1.6165728) q[1];
rz(0.98967057) q[2];
sx q[2];
rz(-2.4151617) q[2];
sx q[2];
rz(-2.5834609) q[2];
rz(0.67047337) q[3];
sx q[3];
rz(-2.5994876) q[3];
sx q[3];
rz(2.2710298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
