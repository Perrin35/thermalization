OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(3.9711877) q[0];
sx q[0];
rz(9.2708099) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(0.33831236) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4727288) q[0];
sx q[0];
rz(-0.40574408) q[0];
sx q[0];
rz(0.91234447) q[0];
rz(-pi) q[1];
rz(-0.36891919) q[2];
sx q[2];
rz(-0.92637617) q[2];
sx q[2];
rz(1.8298139) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.675128) q[1];
sx q[1];
rz(-2.8505241) q[1];
sx q[1];
rz(-0.18578271) q[1];
rz(1.5428513) q[3];
sx q[3];
rz(-2.6290647) q[3];
sx q[3];
rz(2.724109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.14264318) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(-1.1738698) q[2];
rz(-0.075803444) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3409815) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(3.0766292) q[0];
rz(0.57463542) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(-1.8992791) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7466465) q[0];
sx q[0];
rz(-1.6992237) q[0];
sx q[0];
rz(-2.8917679) q[0];
rz(-pi) q[1];
rz(-2.20349) q[2];
sx q[2];
rz(-1.6993076) q[2];
sx q[2];
rz(-3.1046257) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6537135) q[1];
sx q[1];
rz(-2.5839845) q[1];
sx q[1];
rz(0.30456581) q[1];
rz(-pi) q[2];
rz(0.55450704) q[3];
sx q[3];
rz(-2.3508361) q[3];
sx q[3];
rz(-2.8046372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80766455) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(0.57404533) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42049256) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(0.72845355) q[0];
rz(1.4942253) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(2.1247991) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4749878) q[0];
sx q[0];
rz(-3.0525065) q[0];
sx q[0];
rz(2.7276917) q[0];
x q[1];
rz(2.617308) q[2];
sx q[2];
rz(-1.1384083) q[2];
sx q[2];
rz(-2.1081032) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96117561) q[1];
sx q[1];
rz(-0.75913402) q[1];
sx q[1];
rz(1.7708066) q[1];
rz(-pi) q[2];
rz(0.049116491) q[3];
sx q[3];
rz(-1.8521063) q[3];
sx q[3];
rz(1.0055055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3399405) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(2.0920848) q[2];
rz(0.63878757) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(-1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.3291572) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(-1.4720434) q[0];
rz(-2.4064348) q[1];
sx q[1];
rz(-2.3627294) q[1];
sx q[1];
rz(-2.8947815) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8762159) q[0];
sx q[0];
rz(-2.2582158) q[0];
sx q[0];
rz(0.19296293) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6631782) q[2];
sx q[2];
rz(-1.7583499) q[2];
sx q[2];
rz(-2.0699376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88627316) q[1];
sx q[1];
rz(-1.5181932) q[1];
sx q[1];
rz(0.37087755) q[1];
rz(-1.1770583) q[3];
sx q[3];
rz(-0.65440946) q[3];
sx q[3];
rz(2.8670093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4776769) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(1.6003312) q[2];
rz(-2.4345496) q[3];
sx q[3];
rz(-2.0440846) q[3];
sx q[3];
rz(-2.2533806) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8108869) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(1.8141618) q[0];
rz(-1.5785626) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(0.24838233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7029593) q[0];
sx q[0];
rz(-0.54134936) q[0];
sx q[0];
rz(-0.066141733) q[0];
rz(2.9551198) q[2];
sx q[2];
rz(-1.7261793) q[2];
sx q[2];
rz(-1.1764256) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69265428) q[1];
sx q[1];
rz(-0.9874953) q[1];
sx q[1];
rz(1.7765691) q[1];
rz(-2.4928717) q[3];
sx q[3];
rz(-1.3789163) q[3];
sx q[3];
rz(-2.3313525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43626943) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-0.79745897) q[2];
rz(2.752839) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1335063) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(1.3457993) q[0];
rz(-2.0603518) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(0.1246917) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42449441) q[0];
sx q[0];
rz(-1.6569123) q[0];
sx q[0];
rz(2.9647102) q[0];
rz(-pi) q[1];
rz(-0.40839809) q[2];
sx q[2];
rz(-2.4917267) q[2];
sx q[2];
rz(1.0330531) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71619892) q[1];
sx q[1];
rz(-2.4024706) q[1];
sx q[1];
rz(-0.083934099) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3966339) q[3];
sx q[3];
rz(-0.31674851) q[3];
sx q[3];
rz(-1.6572286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6283915) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(-2.52264) q[2];
rz(1.0533054) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(-0.9296023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5597647) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(-1.0429617) q[0];
rz(2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(2.0708864) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571072) q[0];
sx q[0];
rz(-2.6575408) q[0];
sx q[0];
rz(0.0012782106) q[0];
x q[1];
rz(-1.4378909) q[2];
sx q[2];
rz(-1.9153567) q[2];
sx q[2];
rz(-0.29495707) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0106525) q[1];
sx q[1];
rz(-2.3619235) q[1];
sx q[1];
rz(-0.14203771) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0827818) q[3];
sx q[3];
rz(-0.33274129) q[3];
sx q[3];
rz(0.22418338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.36859194) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(0.32361844) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535646) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(-1.2127097) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(0.94747296) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038717) q[0];
sx q[0];
rz(-0.54715711) q[0];
sx q[0];
rz(-1.7659811) q[0];
rz(-pi) q[1];
rz(-2.410789) q[2];
sx q[2];
rz(-2.9060504) q[2];
sx q[2];
rz(1.8159602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.84711134) q[1];
sx q[1];
rz(-1.0080907) q[1];
sx q[1];
rz(-2.860445) q[1];
rz(-pi) q[2];
rz(-0.8074699) q[3];
sx q[3];
rz(-1.8261357) q[3];
sx q[3];
rz(-2.1047999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.56132135) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(-0.46978152) q[2];
rz(-1.3011159) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25205055) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(1.8126194) q[0];
rz(-2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-0.20283094) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348939) q[0];
sx q[0];
rz(-1.5447504) q[0];
sx q[0];
rz(3.1025725) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.780026) q[2];
sx q[2];
rz(-1.3522569) q[2];
sx q[2];
rz(-2.2184559) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.087346615) q[1];
sx q[1];
rz(-3.0294703) q[1];
sx q[1];
rz(2.0263158) q[1];
x q[2];
rz(2.6796954) q[3];
sx q[3];
rz(-1.7089205) q[3];
sx q[3];
rz(-1.671333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41708502) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(0.99651304) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(2.3341808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080169454) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(-2.9598575) q[0];
rz(0.043047992) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(0.28082401) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94175324) q[0];
sx q[0];
rz(-2.0240677) q[0];
sx q[0];
rz(1.0928632) q[0];
rz(-pi) q[1];
rz(-0.21364613) q[2];
sx q[2];
rz(-1.3823969) q[2];
sx q[2];
rz(1.633916) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3481969) q[1];
sx q[1];
rz(-1.5131725) q[1];
sx q[1];
rz(-3.1053931) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0481846) q[3];
sx q[3];
rz(-2.7624353) q[3];
sx q[3];
rz(2.2446333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8250371) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(2.5184856) q[2];
rz(-2.1394219) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(-0.56308693) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6476718) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(2.2676246) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-2.5429824) q[2];
sx q[2];
rz(-2.241588) q[2];
sx q[2];
rz(-0.66551756) q[2];
rz(-0.35691805) q[3];
sx q[3];
rz(-1.1834984) q[3];
sx q[3];
rz(-0.055565861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
