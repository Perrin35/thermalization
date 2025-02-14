OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.35535204) q[0];
sx q[0];
rz(-0.30682895) q[0];
sx q[0];
rz(0.29875779) q[0];
rz(-0.6660676) q[1];
sx q[1];
rz(2.5112285) q[1];
sx q[1];
rz(11.093333) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0399892) q[0];
sx q[0];
rz(-2.371006) q[0];
sx q[0];
rz(-1.5250852) q[0];
x q[1];
rz(3.1236727) q[2];
sx q[2];
rz(-1.8715053) q[2];
sx q[2];
rz(1.199388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.70570488) q[1];
sx q[1];
rz(-1.0748362) q[1];
sx q[1];
rz(0.75261799) q[1];
rz(-pi) q[2];
rz(0.27582081) q[3];
sx q[3];
rz(-1.5032056) q[3];
sx q[3];
rz(1.7662314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1851958) q[2];
sx q[2];
rz(-2.7608725) q[2];
sx q[2];
rz(-1.3107276) q[2];
rz(3.0500566) q[3];
sx q[3];
rz(-2.5444578) q[3];
sx q[3];
rz(-0.53102791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0700584) q[0];
sx q[0];
rz(-1.0907084) q[0];
sx q[0];
rz(-2.6614905) q[0];
rz(1.1473038) q[1];
sx q[1];
rz(-2.5105748) q[1];
sx q[1];
rz(-2.4400585) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6433367) q[0];
sx q[0];
rz(-1.0717046) q[0];
sx q[0];
rz(2.3168111) q[0];
rz(-pi) q[1];
rz(-1.105179) q[2];
sx q[2];
rz(-1.1934416) q[2];
sx q[2];
rz(-1.6060262) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5614723) q[1];
sx q[1];
rz(-1.1271203) q[1];
sx q[1];
rz(0.89718141) q[1];
rz(-0.63752709) q[3];
sx q[3];
rz(-1.3204323) q[3];
sx q[3];
rz(1.7660146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39701617) q[2];
sx q[2];
rz(-1.959356) q[2];
sx q[2];
rz(1.5783295) q[2];
rz(-2.8619316) q[3];
sx q[3];
rz(-0.71664387) q[3];
sx q[3];
rz(2.7320812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.3438943) q[0];
sx q[0];
rz(-2.7293623) q[0];
sx q[0];
rz(-1.4233587) q[0];
rz(0.1238981) q[1];
sx q[1];
rz(-2.2975477) q[1];
sx q[1];
rz(-1.557225) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6829872) q[0];
sx q[0];
rz(-1.2747972) q[0];
sx q[0];
rz(0.83403765) q[0];
x q[1];
rz(-2.2940647) q[2];
sx q[2];
rz(-1.28994) q[2];
sx q[2];
rz(-2.6658863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5277943) q[1];
sx q[1];
rz(-1.8168194) q[1];
sx q[1];
rz(0.08155827) q[1];
x q[2];
rz(-2.2267659) q[3];
sx q[3];
rz(-2.1644691) q[3];
sx q[3];
rz(2.3542924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.065217) q[2];
sx q[2];
rz(-2.4880444) q[2];
sx q[2];
rz(0.93488133) q[2];
rz(-2.8377418) q[3];
sx q[3];
rz(-0.6128208) q[3];
sx q[3];
rz(2.2936308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3728751) q[0];
sx q[0];
rz(-0.77943742) q[0];
sx q[0];
rz(0.26707643) q[0];
rz(0.13485394) q[1];
sx q[1];
rz(-0.57661533) q[1];
sx q[1];
rz(-0.83736247) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98146472) q[0];
sx q[0];
rz(-0.61580393) q[0];
sx q[0];
rz(-1.0976904) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83149399) q[2];
sx q[2];
rz(-2.0250426) q[2];
sx q[2];
rz(-1.9902094) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.048031295) q[1];
sx q[1];
rz(-1.8037027) q[1];
sx q[1];
rz(-1.770411) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6322953) q[3];
sx q[3];
rz(-1.617518) q[3];
sx q[3];
rz(0.72544569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5460633) q[2];
sx q[2];
rz(-0.17543051) q[2];
sx q[2];
rz(2.9626633) q[2];
rz(-1.1692125) q[3];
sx q[3];
rz(-1.7763276) q[3];
sx q[3];
rz(0.21807142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021407481) q[0];
sx q[0];
rz(-2.6299801) q[0];
sx q[0];
rz(0.87274337) q[0];
rz(-1.9841638) q[1];
sx q[1];
rz(-0.87481421) q[1];
sx q[1];
rz(3.1246368) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8411422) q[0];
sx q[0];
rz(-2.4068641) q[0];
sx q[0];
rz(3.1233913) q[0];
rz(2.5418607) q[2];
sx q[2];
rz(-1.5790734) q[2];
sx q[2];
rz(-0.042009609) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8230205) q[1];
sx q[1];
rz(-1.2766395) q[1];
sx q[1];
rz(-2.3579954) q[1];
rz(0.26930974) q[3];
sx q[3];
rz(-1.8519173) q[3];
sx q[3];
rz(0.65015974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.056328) q[2];
sx q[2];
rz(-1.7076098) q[2];
sx q[2];
rz(-2.8223574) q[2];
rz(-0.6790092) q[3];
sx q[3];
rz(-1.7947936) q[3];
sx q[3];
rz(2.3398248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7769258) q[0];
sx q[0];
rz(-0.32061446) q[0];
sx q[0];
rz(-0.6426245) q[0];
rz(1.7721666) q[1];
sx q[1];
rz(-0.6232999) q[1];
sx q[1];
rz(-2.1646037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1851121) q[0];
sx q[0];
rz(-1.4727043) q[0];
sx q[0];
rz(0.12382671) q[0];
rz(1.0426635) q[2];
sx q[2];
rz(-2.0181877) q[2];
sx q[2];
rz(-1.7698947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.381666) q[1];
sx q[1];
rz(-1.5846415) q[1];
sx q[1];
rz(-2.2626876) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6518408) q[3];
sx q[3];
rz(-0.5872927) q[3];
sx q[3];
rz(1.0451661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52648181) q[2];
sx q[2];
rz(-1.4075764) q[2];
sx q[2];
rz(1.2283481) q[2];
rz(2.2743716) q[3];
sx q[3];
rz(-0.4129748) q[3];
sx q[3];
rz(0.76798463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7277471) q[0];
sx q[0];
rz(-0.21738805) q[0];
sx q[0];
rz(-0.17284285) q[0];
rz(2.3364283) q[1];
sx q[1];
rz(-2.5925437) q[1];
sx q[1];
rz(3.0056312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7522536) q[0];
sx q[0];
rz(-0.87119188) q[0];
sx q[0];
rz(-2.9851578) q[0];
rz(1.4589086) q[2];
sx q[2];
rz(-0.75190836) q[2];
sx q[2];
rz(-2.65154) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0522642) q[1];
sx q[1];
rz(-1.5517527) q[1];
sx q[1];
rz(0.11472265) q[1];
rz(-pi) q[2];
rz(0.32973955) q[3];
sx q[3];
rz(-1.6981594) q[3];
sx q[3];
rz(-1.691526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.92676306) q[2];
sx q[2];
rz(-1.5952933) q[2];
sx q[2];
rz(2.9570441) q[2];
rz(-0.18937011) q[3];
sx q[3];
rz(-0.53818494) q[3];
sx q[3];
rz(2.2034933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9823343) q[0];
sx q[0];
rz(-2.8026447) q[0];
sx q[0];
rz(0.92639297) q[0];
rz(0.039208086) q[1];
sx q[1];
rz(-2.4640633) q[1];
sx q[1];
rz(1.0367941) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2873357) q[0];
sx q[0];
rz(-1.2764036) q[0];
sx q[0];
rz(-2.9578985) q[0];
rz(-2.7095058) q[2];
sx q[2];
rz(-1.7361904) q[2];
sx q[2];
rz(1.337932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7651789) q[1];
sx q[1];
rz(-2.661099) q[1];
sx q[1];
rz(-2.5960467) q[1];
rz(2.751334) q[3];
sx q[3];
rz(-1.3969799) q[3];
sx q[3];
rz(2.9040608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2490273) q[2];
sx q[2];
rz(-1.1250863) q[2];
sx q[2];
rz(-0.79891515) q[2];
rz(-0.51698452) q[3];
sx q[3];
rz(-0.40701443) q[3];
sx q[3];
rz(-2.5911205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35140458) q[0];
sx q[0];
rz(-3.1365972) q[0];
sx q[0];
rz(1.5193526) q[0];
rz(1.5478569) q[1];
sx q[1];
rz(-2.3212815) q[1];
sx q[1];
rz(2.4511852) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6192157) q[0];
sx q[0];
rz(-1.7937614) q[0];
sx q[0];
rz(0.35943835) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79924714) q[2];
sx q[2];
rz(-1.4991781) q[2];
sx q[2];
rz(-1.3992753) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9195329) q[1];
sx q[1];
rz(-1.0708978) q[1];
sx q[1];
rz(-0.57603194) q[1];
rz(-2.1097095) q[3];
sx q[3];
rz(-1.6953985) q[3];
sx q[3];
rz(-2.3923158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4385628) q[2];
sx q[2];
rz(-2.5327693) q[2];
sx q[2];
rz(-2.1717066) q[2];
rz(0.25740933) q[3];
sx q[3];
rz(-2.9195547) q[3];
sx q[3];
rz(1.359587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6064706) q[0];
sx q[0];
rz(-0.73000014) q[0];
sx q[0];
rz(3.0560793) q[0];
rz(-2.8657148) q[1];
sx q[1];
rz(-0.62092263) q[1];
sx q[1];
rz(1.9524908) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082423) q[0];
sx q[0];
rz(-1.9353284) q[0];
sx q[0];
rz(-0.43594283) q[0];
rz(-pi) q[1];
rz(0.027935426) q[2];
sx q[2];
rz(-0.31012529) q[2];
sx q[2];
rz(-1.5903697) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.72451776) q[1];
sx q[1];
rz(-2.7226884) q[1];
sx q[1];
rz(-1.8663097) q[1];
x q[2];
rz(2.7932634) q[3];
sx q[3];
rz(-1.4831721) q[3];
sx q[3];
rz(1.5867811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4756061) q[2];
sx q[2];
rz(-0.87916547) q[2];
sx q[2];
rz(2.6574668) q[2];
rz(0.13949805) q[3];
sx q[3];
rz(-0.77553427) q[3];
sx q[3];
rz(2.9805396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6763247) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(1.4576661) q[1];
sx q[1];
rz(-1.6438345) q[1];
sx q[1];
rz(1.8297292) q[1];
rz(-2.8756782) q[2];
sx q[2];
rz(-2.4752046) q[2];
sx q[2];
rz(2.2023329) q[2];
rz(2.3845354) q[3];
sx q[3];
rz(-1.591605) q[3];
sx q[3];
rz(0.021416728) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
