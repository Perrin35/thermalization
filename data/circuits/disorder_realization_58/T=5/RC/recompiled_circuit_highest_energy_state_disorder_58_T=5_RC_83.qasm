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
rz(0.27920029) q[0];
sx q[0];
rz(-1.2503799) q[0];
sx q[0];
rz(3.1412636) q[0];
rz(-2.012913) q[1];
sx q[1];
rz(-2.0560775) q[1];
sx q[1];
rz(0.20927277) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64855382) q[0];
sx q[0];
rz(-1.4161515) q[0];
sx q[0];
rz(-1.3376608) q[0];
rz(0.26165719) q[2];
sx q[2];
rz(-2.5664866) q[2];
sx q[2];
rz(1.9272137) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55587253) q[1];
sx q[1];
rz(-1.861708) q[1];
sx q[1];
rz(-2.8147999) q[1];
rz(-pi) q[2];
rz(-1.3253778) q[3];
sx q[3];
rz(-1.8031464) q[3];
sx q[3];
rz(1.7690675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2613232) q[2];
sx q[2];
rz(-1.3618733) q[2];
sx q[2];
rz(-0.23639354) q[2];
rz(-2.7586625) q[3];
sx q[3];
rz(-2.0536486) q[3];
sx q[3];
rz(-1.6898164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3799389) q[0];
sx q[0];
rz(-2.8528657) q[0];
sx q[0];
rz(2.5987103) q[0];
rz(-2.701243) q[1];
sx q[1];
rz(-1.1232702) q[1];
sx q[1];
rz(2.0139258) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714345) q[0];
sx q[0];
rz(-1.9160524) q[0];
sx q[0];
rz(-2.4717114) q[0];
x q[1];
rz(1.3593349) q[2];
sx q[2];
rz(-0.64675719) q[2];
sx q[2];
rz(0.46403592) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.013036916) q[1];
sx q[1];
rz(-2.3686976) q[1];
sx q[1];
rz(-1.3580639) q[1];
rz(2.4324969) q[3];
sx q[3];
rz(-1.2768687) q[3];
sx q[3];
rz(-0.16277587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.27309624) q[2];
sx q[2];
rz(-2.345583) q[2];
sx q[2];
rz(-1.8491171) q[2];
rz(0.83186197) q[3];
sx q[3];
rz(-1.4757194) q[3];
sx q[3];
rz(2.7003435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.454994) q[0];
sx q[0];
rz(-0.43231493) q[0];
sx q[0];
rz(-1.5469714) q[0];
rz(3.1349685) q[1];
sx q[1];
rz(-0.5358271) q[1];
sx q[1];
rz(0.63820249) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4362559) q[0];
sx q[0];
rz(-0.47644189) q[0];
sx q[0];
rz(1.492155) q[0];
rz(-0.77000051) q[2];
sx q[2];
rz(-1.8714617) q[2];
sx q[2];
rz(-1.2141924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6276826) q[1];
sx q[1];
rz(-1.9238775) q[1];
sx q[1];
rz(-0.38796918) q[1];
rz(-pi) q[2];
rz(1.9157103) q[3];
sx q[3];
rz(-2.5673491) q[3];
sx q[3];
rz(1.739603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39614761) q[2];
sx q[2];
rz(-0.026693176) q[2];
sx q[2];
rz(0.8566345) q[2];
rz(0.73532295) q[3];
sx q[3];
rz(-1.4152799) q[3];
sx q[3];
rz(2.0935238) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.624144) q[0];
sx q[0];
rz(-0.051001661) q[0];
sx q[0];
rz(0.93358246) q[0];
rz(2.3816662) q[1];
sx q[1];
rz(-2.3257207) q[1];
sx q[1];
rz(1.5911721) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5869414) q[0];
sx q[0];
rz(-1.6593895) q[0];
sx q[0];
rz(2.8604126) q[0];
rz(0.43689219) q[2];
sx q[2];
rz(-1.6342739) q[2];
sx q[2];
rz(-1.2530099) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62867579) q[1];
sx q[1];
rz(-2.4077291) q[1];
sx q[1];
rz(0.34725125) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.891734) q[3];
sx q[3];
rz(-1.7267531) q[3];
sx q[3];
rz(0.74710959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.4178702) q[2];
sx q[2];
rz(-1.7524065) q[2];
sx q[2];
rz(2.1088481) q[2];
rz(-0.51747733) q[3];
sx q[3];
rz(-1.7101219) q[3];
sx q[3];
rz(-1.8364068) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0378157) q[0];
sx q[0];
rz(-3.0396437) q[0];
sx q[0];
rz(1.284449) q[0];
rz(2.2880554) q[1];
sx q[1];
rz(-1.6910911) q[1];
sx q[1];
rz(2.6654713) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0697107) q[0];
sx q[0];
rz(-0.5997385) q[0];
sx q[0];
rz(1.016933) q[0];
x q[1];
rz(-2.7760162) q[2];
sx q[2];
rz(-1.3028868) q[2];
sx q[2];
rz(0.6838201) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13083982) q[1];
sx q[1];
rz(-0.71341842) q[1];
sx q[1];
rz(-2.6793007) q[1];
rz(1.4912603) q[3];
sx q[3];
rz(-2.400268) q[3];
sx q[3];
rz(2.0052412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7242754) q[2];
sx q[2];
rz(-2.3099835) q[2];
sx q[2];
rz(-1.6012491) q[2];
rz(0.075751461) q[3];
sx q[3];
rz(-1.3906393) q[3];
sx q[3];
rz(0.72995228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1002355) q[0];
sx q[0];
rz(-0.55158177) q[0];
sx q[0];
rz(1.7919737) q[0];
rz(-1.8655818) q[1];
sx q[1];
rz(-0.7834692) q[1];
sx q[1];
rz(-1.0412201) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2305377) q[0];
sx q[0];
rz(-1.7887702) q[0];
sx q[0];
rz(2.387052) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5346933) q[2];
sx q[2];
rz(-2.3182333) q[2];
sx q[2];
rz(-0.90561101) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1043084) q[1];
sx q[1];
rz(-1.3194024) q[1];
sx q[1];
rz(-2.7576552) q[1];
rz(-pi) q[2];
rz(-1.8550664) q[3];
sx q[3];
rz(-0.83300029) q[3];
sx q[3];
rz(-2.7018911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2714046) q[2];
sx q[2];
rz(-1.1689593) q[2];
sx q[2];
rz(2.1427587) q[2];
rz(-2.6073604) q[3];
sx q[3];
rz(-1.1778919) q[3];
sx q[3];
rz(-1.6629201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44487318) q[0];
sx q[0];
rz(-2.1598926) q[0];
sx q[0];
rz(-1.5848507) q[0];
rz(1.0818256) q[1];
sx q[1];
rz(-1.6422681) q[1];
sx q[1];
rz(-0.36344847) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9173319) q[0];
sx q[0];
rz(-1.5646569) q[0];
sx q[0];
rz(0.21256769) q[0];
rz(-pi) q[1];
rz(2.8055023) q[2];
sx q[2];
rz(-0.31340965) q[2];
sx q[2];
rz(1.6801419) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9103437) q[1];
sx q[1];
rz(-1.1414613) q[1];
sx q[1];
rz(-2.3554159) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0814675) q[3];
sx q[3];
rz(-0.31144588) q[3];
sx q[3];
rz(-1.2490727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58793679) q[2];
sx q[2];
rz(-2.9755972) q[2];
sx q[2];
rz(2.0641649) q[2];
rz(-1.59498) q[3];
sx q[3];
rz(-1.8190106) q[3];
sx q[3];
rz(-0.67302978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3308291) q[0];
sx q[0];
rz(-1.6442693) q[0];
sx q[0];
rz(-0.28054917) q[0];
rz(0.55889091) q[1];
sx q[1];
rz(-2.2160896) q[1];
sx q[1];
rz(-1.2116609) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9292363) q[0];
sx q[0];
rz(-0.87221891) q[0];
sx q[0];
rz(0.86232604) q[0];
rz(-pi) q[1];
rz(-2.6490491) q[2];
sx q[2];
rz(-1.768868) q[2];
sx q[2];
rz(-2.1576674) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9338464) q[1];
sx q[1];
rz(-2.3504815) q[1];
sx q[1];
rz(1.6027811) q[1];
rz(1.7688076) q[3];
sx q[3];
rz(-1.5757547) q[3];
sx q[3];
rz(1.229242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9676548) q[2];
sx q[2];
rz(-1.8123947) q[2];
sx q[2];
rz(0.8052513) q[2];
rz(2.0179181) q[3];
sx q[3];
rz(-1.1403964) q[3];
sx q[3];
rz(-2.8203188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6174018) q[0];
sx q[0];
rz(-0.62365714) q[0];
sx q[0];
rz(-2.6362841) q[0];
rz(2.4303923) q[1];
sx q[1];
rz(-2.2732747) q[1];
sx q[1];
rz(2.0288846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8998535) q[0];
sx q[0];
rz(-2.02075) q[0];
sx q[0];
rz(-1.8983272) q[0];
rz(-pi) q[1];
rz(2.4975487) q[2];
sx q[2];
rz(-0.13453787) q[2];
sx q[2];
rz(0.74468871) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0074626) q[1];
sx q[1];
rz(-0.45744823) q[1];
sx q[1];
rz(1.8602954) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5212413) q[3];
sx q[3];
rz(-0.55618868) q[3];
sx q[3];
rz(2.9002528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6342371) q[2];
sx q[2];
rz(-1.2398182) q[2];
sx q[2];
rz(-0.235454) q[2];
rz(2.7638451) q[3];
sx q[3];
rz(-2.7538959) q[3];
sx q[3];
rz(3.0915251) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3560051) q[0];
sx q[0];
rz(-1.8536785) q[0];
sx q[0];
rz(0.054314286) q[0];
rz(2.3771225) q[1];
sx q[1];
rz(-0.46934325) q[1];
sx q[1];
rz(2.4321709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5886836) q[0];
sx q[0];
rz(-1.3257265) q[0];
sx q[0];
rz(-0.10677307) q[0];
rz(-pi) q[1];
rz(1.7058175) q[2];
sx q[2];
rz(-1.8293937) q[2];
sx q[2];
rz(2.1136346) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.874234) q[1];
sx q[1];
rz(-2.5399985) q[1];
sx q[1];
rz(-0.20341058) q[1];
rz(-1.7001496) q[3];
sx q[3];
rz(-1.23063) q[3];
sx q[3];
rz(2.8486855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28598049) q[2];
sx q[2];
rz(-1.6667655) q[2];
sx q[2];
rz(-0.093185514) q[2];
rz(-2.5423673) q[3];
sx q[3];
rz(-0.56120509) q[3];
sx q[3];
rz(-0.78527251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7080606) q[0];
sx q[0];
rz(-2.1716433) q[0];
sx q[0];
rz(-0.94950983) q[0];
rz(3.0769729) q[1];
sx q[1];
rz(-3.0977991) q[1];
sx q[1];
rz(-2.480712) q[1];
rz(-1.6914037) q[2];
sx q[2];
rz(-2.5868138) q[2];
sx q[2];
rz(-0.13704337) q[2];
rz(0.6653847) q[3];
sx q[3];
rz(-1.7177842) q[3];
sx q[3];
rz(-0.73242889) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
