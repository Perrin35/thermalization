OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(-1.7927875) q[1];
sx q[1];
rz(-0.92372149) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38102725) q[0];
sx q[0];
rz(-1.1228704) q[0];
sx q[0];
rz(0.3737803) q[0];
rz(0.14416868) q[2];
sx q[2];
rz(-1.2906133) q[2];
sx q[2];
rz(-0.35137128) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10780653) q[1];
sx q[1];
rz(-2.4596655) q[1];
sx q[1];
rz(-2.4297907) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5039623) q[3];
sx q[3];
rz(-2.5507567) q[3];
sx q[3];
rz(1.5096111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2279921) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(2.9795734) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.0733923) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.4555567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14917063) q[0];
sx q[0];
rz(-2.143321) q[0];
sx q[0];
rz(0.19897977) q[0];
x q[1];
rz(-2.7669719) q[2];
sx q[2];
rz(-1.4943559) q[2];
sx q[2];
rz(2.3945216) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96869722) q[1];
sx q[1];
rz(-0.52297938) q[1];
sx q[1];
rz(1.5437267) q[1];
x q[2];
rz(-1.5973813) q[3];
sx q[3];
rz(-1.9823325) q[3];
sx q[3];
rz(2.6708024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.3519752) q[2];
rz(-0.18243608) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(-2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75333726) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(-1.9690537) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81329936) q[0];
sx q[0];
rz(-1.0070224) q[0];
sx q[0];
rz(1.9378807) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8981947) q[2];
sx q[2];
rz(-0.84257579) q[2];
sx q[2];
rz(-0.74795216) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.064627083) q[1];
sx q[1];
rz(-0.58276999) q[1];
sx q[1];
rz(-1.1175734) q[1];
rz(-pi) q[2];
rz(-1.8215239) q[3];
sx q[3];
rz(-1.6276976) q[3];
sx q[3];
rz(-2.1123321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(0.91119901) q[2];
rz(-0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.38055414) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(0.52350837) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9317755) q[0];
sx q[0];
rz(-2.1486001) q[0];
sx q[0];
rz(1.125976) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8977215) q[2];
sx q[2];
rz(-0.3393617) q[2];
sx q[2];
rz(1.3432168) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7586786) q[1];
sx q[1];
rz(-2.3944693) q[1];
sx q[1];
rz(-1.1927356) q[1];
x q[2];
rz(1.773049) q[3];
sx q[3];
rz(-2.5407255) q[3];
sx q[3];
rz(-0.7522538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48150912) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(-2.2619757) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-0.99194828) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.49682) q[0];
sx q[0];
rz(-0.20475514) q[0];
sx q[0];
rz(-2.5909008) q[0];
x q[1];
rz(1.5993824) q[2];
sx q[2];
rz(-2.8386142) q[2];
sx q[2];
rz(-2.6945393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3494898) q[1];
sx q[1];
rz(-1.2621242) q[1];
sx q[1];
rz(-1.5285138) q[1];
rz(-1.7080073) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(-1.8828132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(0.27080718) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4145684) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(-1.4250925) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3544918) q[0];
sx q[0];
rz(-1.0757425) q[0];
sx q[0];
rz(1.8255193) q[0];
rz(-1.9082597) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(-2.595682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.35364321) q[1];
sx q[1];
rz(-2.0197581) q[1];
sx q[1];
rz(0.073972703) q[1];
rz(-2.1720042) q[3];
sx q[3];
rz(-2.0242656) q[3];
sx q[3];
rz(-0.78905247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(-2.5777204) q[2];
rz(-0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5086223) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(2.7222743) q[0];
rz(-1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(2.3197876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069698378) q[0];
sx q[0];
rz(-2.0928203) q[0];
sx q[0];
rz(-0.72899039) q[0];
rz(-pi) q[1];
rz(2.6452438) q[2];
sx q[2];
rz(-0.90663547) q[2];
sx q[2];
rz(-1.6830483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0710443) q[1];
sx q[1];
rz(-0.76862915) q[1];
sx q[1];
rz(-2.8119836) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5542198) q[3];
sx q[3];
rz(-1.2578739) q[3];
sx q[3];
rz(-1.3164933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3372779) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(-0.129536) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(2.3186671) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(1.3051422) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6179498) q[0];
sx q[0];
rz(-2.4542913) q[0];
sx q[0];
rz(0.53074093) q[0];
rz(-pi) q[1];
rz(0.70456409) q[2];
sx q[2];
rz(-1.8578055) q[2];
sx q[2];
rz(-1.9412083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23755632) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(2.0380286) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2315024) q[3];
sx q[3];
rz(-1.5486071) q[3];
sx q[3];
rz(2.0458178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(-1.0772609) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(2.819678) q[0];
rz(1.6053258) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-0.70294356) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020333175) q[0];
sx q[0];
rz(-2.5998305) q[0];
sx q[0];
rz(2.3938177) q[0];
x q[1];
rz(-3.1019194) q[2];
sx q[2];
rz(-2.0797605) q[2];
sx q[2];
rz(-0.85658011) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1785537) q[1];
sx q[1];
rz(-2.5759765) q[1];
sx q[1];
rz(-1.8815243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1144936) q[3];
sx q[3];
rz(-2.4745686) q[3];
sx q[3];
rz(-1.4240571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(-0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777578) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(0.46863619) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2221453) q[0];
sx q[0];
rz(-1.7788017) q[0];
sx q[0];
rz(-0.39214765) q[0];
x q[1];
rz(-0.87601985) q[2];
sx q[2];
rz(-0.48601905) q[2];
sx q[2];
rz(1.8234058) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.66099) q[1];
sx q[1];
rz(-2.8669679) q[1];
sx q[1];
rz(1.8250699) q[1];
x q[2];
rz(2.7077984) q[3];
sx q[3];
rz(-1.9926096) q[3];
sx q[3];
rz(-2.1123561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.1432077) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951915) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(2.519683) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(2.2767699) q[2];
sx q[2];
rz(-2.13158) q[2];
sx q[2];
rz(-2.0958015) q[2];
rz(-3.0631089) q[3];
sx q[3];
rz(-2.2208636) q[3];
sx q[3];
rz(0.29490864) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
