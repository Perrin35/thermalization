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
rz(2.3693585) q[0];
sx q[0];
rz(-1.9263664) q[0];
sx q[0];
rz(-1.9908494) q[0];
rz(2.8590705) q[1];
sx q[1];
rz(-1.1257659) q[1];
sx q[1];
rz(-1.1123302) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.756506) q[0];
sx q[0];
rz(-1.6876564) q[0];
sx q[0];
rz(-2.8782842) q[0];
rz(-pi) q[1];
rz(-1.7409647) q[2];
sx q[2];
rz(-2.552816) q[2];
sx q[2];
rz(2.6937304) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9196786) q[1];
sx q[1];
rz(-2.3828607) q[1];
sx q[1];
rz(-0.83617156) q[1];
rz(2.2745773) q[3];
sx q[3];
rz(-0.91546446) q[3];
sx q[3];
rz(1.0176147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7925966) q[2];
sx q[2];
rz(-1.0745606) q[2];
sx q[2];
rz(1.7462771) q[2];
rz(-3.0830749) q[3];
sx q[3];
rz(-1.7250215) q[3];
sx q[3];
rz(1.831656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42404744) q[0];
sx q[0];
rz(-0.40322867) q[0];
sx q[0];
rz(2.6918217) q[0];
rz(0.46478081) q[1];
sx q[1];
rz(-1.8315146) q[1];
sx q[1];
rz(-0.25258499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3656986) q[0];
sx q[0];
rz(-1.5868385) q[0];
sx q[0];
rz(-1.2195862) q[0];
rz(-pi) q[1];
rz(2.4812883) q[2];
sx q[2];
rz(-1.2129158) q[2];
sx q[2];
rz(1.7972657) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8014268) q[1];
sx q[1];
rz(-0.28160843) q[1];
sx q[1];
rz(-1.5605218) q[1];
x q[2];
rz(-0.95460598) q[3];
sx q[3];
rz(-2.2807215) q[3];
sx q[3];
rz(1.4675642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2437336) q[2];
sx q[2];
rz(-0.89577883) q[2];
sx q[2];
rz(-0.016156999) q[2];
rz(0.92205087) q[3];
sx q[3];
rz(-2.1478896) q[3];
sx q[3];
rz(3.0097358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048412662) q[0];
sx q[0];
rz(-1.6572297) q[0];
sx q[0];
rz(-0.68159252) q[0];
rz(0.59773481) q[1];
sx q[1];
rz(-1.6860516) q[1];
sx q[1];
rz(2.8173503) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77406835) q[0];
sx q[0];
rz(-1.9107959) q[0];
sx q[0];
rz(2.6794479) q[0];
rz(-pi) q[1];
rz(-0.50679548) q[2];
sx q[2];
rz(-1.629771) q[2];
sx q[2];
rz(0.61117327) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8719667) q[1];
sx q[1];
rz(-1.5287279) q[1];
sx q[1];
rz(1.4021238) q[1];
rz(-pi) q[2];
rz(1.3310166) q[3];
sx q[3];
rz(-0.62285715) q[3];
sx q[3];
rz(1.1774802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9165667) q[2];
sx q[2];
rz(-2.7498507) q[2];
sx q[2];
rz(-2.0908835) q[2];
rz(2.4010036) q[3];
sx q[3];
rz(-1.8480443) q[3];
sx q[3];
rz(-0.061323969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4968313) q[0];
sx q[0];
rz(-0.83503857) q[0];
sx q[0];
rz(-2.187425) q[0];
rz(1.1369811) q[1];
sx q[1];
rz(-1.6258207) q[1];
sx q[1];
rz(-0.7581996) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6369259) q[0];
sx q[0];
rz(-1.6109137) q[0];
sx q[0];
rz(1.841943) q[0];
rz(-1.3761282) q[2];
sx q[2];
rz(-0.83651453) q[2];
sx q[2];
rz(2.0912974) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7247347) q[1];
sx q[1];
rz(-1.4730037) q[1];
sx q[1];
rz(-2.2657402) q[1];
rz(-0.73885804) q[3];
sx q[3];
rz(-1.5631873) q[3];
sx q[3];
rz(1.348902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.053146426) q[2];
sx q[2];
rz(-0.55067486) q[2];
sx q[2];
rz(-1.3147563) q[2];
rz(-2.8407319) q[3];
sx q[3];
rz(-1.8010537) q[3];
sx q[3];
rz(0.68527591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9984197) q[0];
sx q[0];
rz(-1.5851861) q[0];
sx q[0];
rz(1.6402624) q[0];
rz(-2.7550664) q[1];
sx q[1];
rz(-2.0730348) q[1];
sx q[1];
rz(-1.5949257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0536755) q[0];
sx q[0];
rz(-2.2639416) q[0];
sx q[0];
rz(1.7960834) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1793381) q[2];
sx q[2];
rz(-1.3013162) q[2];
sx q[2];
rz(0.4294006) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3749126) q[1];
sx q[1];
rz(-2.0030352) q[1];
sx q[1];
rz(1.47827) q[1];
rz(-3.0647072) q[3];
sx q[3];
rz(-2.5664483) q[3];
sx q[3];
rz(-2.157876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.74653643) q[2];
sx q[2];
rz(-0.48607963) q[2];
sx q[2];
rz(-1.0699832) q[2];
rz(0.88548958) q[3];
sx q[3];
rz(-1.0642137) q[3];
sx q[3];
rz(1.1522256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-2.2162061) q[0];
sx q[0];
rz(-0.33846551) q[0];
sx q[0];
rz(-1.0534519) q[0];
rz(0.18928754) q[1];
sx q[1];
rz(-0.82866755) q[1];
sx q[1];
rz(1.6493571) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5439592) q[0];
sx q[0];
rz(-1.9372371) q[0];
sx q[0];
rz(-0.14800114) q[0];
rz(-0.97215488) q[2];
sx q[2];
rz(-0.36566662) q[2];
sx q[2];
rz(-0.36520457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.511586) q[1];
sx q[1];
rz(-1.6701856) q[1];
sx q[1];
rz(-0.50390999) q[1];
rz(-1.5297278) q[3];
sx q[3];
rz(-0.76309915) q[3];
sx q[3];
rz(1.3535318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3377127) q[2];
sx q[2];
rz(-2.1264075) q[2];
sx q[2];
rz(2.0886776) q[2];
rz(-0.11689154) q[3];
sx q[3];
rz(-2.0785073) q[3];
sx q[3];
rz(-1.6803928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9601124) q[0];
sx q[0];
rz(-0.57992613) q[0];
sx q[0];
rz(0.54452288) q[0];
rz(-0.25360423) q[1];
sx q[1];
rz(-2.8896152) q[1];
sx q[1];
rz(2.8869218) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47557467) q[0];
sx q[0];
rz(-0.072192497) q[0];
sx q[0];
rz(1.6572787) q[0];
x q[1];
rz(2.9355132) q[2];
sx q[2];
rz(-2.7408618) q[2];
sx q[2];
rz(1.3361267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0510301) q[1];
sx q[1];
rz(-1.6335618) q[1];
sx q[1];
rz(-2.323137) q[1];
rz(-2.5768877) q[3];
sx q[3];
rz(-2.2484591) q[3];
sx q[3];
rz(2.0155747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2250681) q[2];
sx q[2];
rz(-2.2038286) q[2];
sx q[2];
rz(2.7226105) q[2];
rz(0.032912832) q[3];
sx q[3];
rz(-1.2337647) q[3];
sx q[3];
rz(1.1121496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58881775) q[0];
sx q[0];
rz(-0.7426312) q[0];
sx q[0];
rz(1.8137929) q[0];
rz(2.1090419) q[1];
sx q[1];
rz(-1.3464709) q[1];
sx q[1];
rz(-0.58951497) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.179993) q[0];
sx q[0];
rz(-2.364453) q[0];
sx q[0];
rz(-0.83013089) q[0];
rz(1.8399182) q[2];
sx q[2];
rz(-1.8523714) q[2];
sx q[2];
rz(0.21965227) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3302916) q[1];
sx q[1];
rz(-1.9988235) q[1];
sx q[1];
rz(2.2045662) q[1];
rz(2.74624) q[3];
sx q[3];
rz(-2.0369923) q[3];
sx q[3];
rz(-0.485093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0716268) q[2];
sx q[2];
rz(-2.4667141) q[2];
sx q[2];
rz(-0.91342941) q[2];
rz(-0.45201284) q[3];
sx q[3];
rz(-1.2830696) q[3];
sx q[3];
rz(-1.7970596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2695059) q[0];
sx q[0];
rz(-2.1156613) q[0];
sx q[0];
rz(-1.5103229) q[0];
rz(-1.3144846) q[1];
sx q[1];
rz(-1.038082) q[1];
sx q[1];
rz(2.7297535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56758948) q[0];
sx q[0];
rz(-1.2943983) q[0];
sx q[0];
rz(0.92692356) q[0];
x q[1];
rz(0.31899231) q[2];
sx q[2];
rz(-1.5249494) q[2];
sx q[2];
rz(0.09825966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50166124) q[1];
sx q[1];
rz(-1.2218214) q[1];
sx q[1];
rz(0.56569143) q[1];
rz(2.4437636) q[3];
sx q[3];
rz(-2.6549454) q[3];
sx q[3];
rz(2.6634288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1371548) q[2];
sx q[2];
rz(-1.4919446) q[2];
sx q[2];
rz(-0.23492661) q[2];
rz(-0.91819417) q[3];
sx q[3];
rz(-2.4166959) q[3];
sx q[3];
rz(-1.0880067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.039577) q[0];
sx q[0];
rz(-0.35440847) q[0];
sx q[0];
rz(1.2603731) q[0];
rz(-1.4541939) q[1];
sx q[1];
rz(-1.4581542) q[1];
sx q[1];
rz(-1.4448602) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8704027) q[0];
sx q[0];
rz(-1.1605485) q[0];
sx q[0];
rz(-2.6543186) q[0];
rz(-pi) q[1];
rz(2.2606662) q[2];
sx q[2];
rz(-1.147715) q[2];
sx q[2];
rz(1.0151052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.48669882) q[1];
sx q[1];
rz(-1.4461262) q[1];
sx q[1];
rz(2.423684) q[1];
rz(-pi) q[2];
rz(-1.904018) q[3];
sx q[3];
rz(-1.929627) q[3];
sx q[3];
rz(-2.8976208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5626278) q[2];
sx q[2];
rz(-0.29890385) q[2];
sx q[2];
rz(-3.01801) q[2];
rz(2.852671) q[3];
sx q[3];
rz(-1.2762504) q[3];
sx q[3];
rz(-2.6821274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73744437) q[0];
sx q[0];
rz(-1.3442208) q[0];
sx q[0];
rz(-1.5303045) q[0];
rz(2.1743446) q[1];
sx q[1];
rz(-0.9728685) q[1];
sx q[1];
rz(-2.2546993) q[1];
rz(-1.6254454) q[2];
sx q[2];
rz(-2.1693632) q[2];
sx q[2];
rz(-2.7827435) q[2];
rz(-0.88913118) q[3];
sx q[3];
rz(-2.0778772) q[3];
sx q[3];
rz(-2.7419013) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
