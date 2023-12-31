OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(3.1728035) q[0];
sx q[0];
rz(6.7682545) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(2.7273942) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.30487) q[0];
sx q[0];
rz(-0.49855907) q[0];
sx q[0];
rz(-2.3303633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0010927) q[2];
sx q[2];
rz(-1.6506519) q[2];
sx q[2];
rz(2.9542838) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84250433) q[1];
sx q[1];
rz(-2.7111004) q[1];
sx q[1];
rz(-1.4166142) q[1];
rz(-pi) q[2];
rz(1.4685417) q[3];
sx q[3];
rz(-1.8342606) q[3];
sx q[3];
rz(1.2276358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2177314) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(0.031575354) q[2];
rz(-1.8850373) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(-2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8388222) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(-0.1698499) q[0];
rz(2.4376712) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-2.6020715) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182506) q[0];
sx q[0];
rz(-1.6174416) q[0];
sx q[0];
rz(-1.1211066) q[0];
rz(-pi) q[1];
rz(-1.1621446) q[2];
sx q[2];
rz(-0.89102972) q[2];
sx q[2];
rz(-1.0811999) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.19903781) q[1];
sx q[1];
rz(-0.53831646) q[1];
sx q[1];
rz(-1.1953137) q[1];
rz(2.2803454) q[3];
sx q[3];
rz(-2.0112787) q[3];
sx q[3];
rz(1.7064106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66723055) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(0.24307069) q[2];
rz(0.66611755) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(-1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3617525) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(0.4483805) q[0];
rz(-1.7547296) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(-0.2562491) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6305144) q[0];
sx q[0];
rz(-0.96195463) q[0];
sx q[0];
rz(-1.8947253) q[0];
rz(-pi) q[1];
rz(-2.3149812) q[2];
sx q[2];
rz(-0.6140784) q[2];
sx q[2];
rz(2.5342864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5445404) q[1];
sx q[1];
rz(-1.2433194) q[1];
sx q[1];
rz(-2.271133) q[1];
rz(-pi) q[2];
rz(0.17685299) q[3];
sx q[3];
rz(-1.5336509) q[3];
sx q[3];
rz(-1.1249441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8024575) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(-1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9280076) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(-1.150594) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(2.4096699) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5314732) q[0];
sx q[0];
rz(-1.1024794) q[0];
sx q[0];
rz(2.5584695) q[0];
rz(-0.133693) q[2];
sx q[2];
rz(-0.72172726) q[2];
sx q[2];
rz(-2.0267817) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5922028) q[1];
sx q[1];
rz(-1.376774) q[1];
sx q[1];
rz(2.8787896) q[1];
rz(1.4452403) q[3];
sx q[3];
rz(-1.7896277) q[3];
sx q[3];
rz(2.9343176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(2.990492) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.2623825) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(0.79777065) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20194963) q[0];
sx q[0];
rz(-1.573283) q[0];
sx q[0];
rz(-1.4506838) q[0];
x q[1];
rz(1.351379) q[2];
sx q[2];
rz(-1.2741538) q[2];
sx q[2];
rz(-0.23362939) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7792015) q[1];
sx q[1];
rz(-2.5070842) q[1];
sx q[1];
rz(1.4125376) q[1];
x q[2];
rz(-2.5161414) q[3];
sx q[3];
rz(-0.52895412) q[3];
sx q[3];
rz(1.4048502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.34565869) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(-0.30203715) q[2];
rz(1.1473514) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(2.5938477) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40925947) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(1.9858032) q[0];
rz(2.0571158) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(-3.0715122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1228186) q[0];
sx q[0];
rz(-2.9111324) q[0];
sx q[0];
rz(-2.3019058) q[0];
x q[1];
rz(2.3577865) q[2];
sx q[2];
rz(-1.98181) q[2];
sx q[2];
rz(-0.54223261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49680432) q[1];
sx q[1];
rz(-1.7323238) q[1];
sx q[1];
rz(0.5715538) q[1];
rz(-pi) q[2];
rz(-1.4962247) q[3];
sx q[3];
rz(-1.5354904) q[3];
sx q[3];
rz(1.1535742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8391116) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(-1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550585) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(0.85987464) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(0.0079356114) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7769653) q[0];
sx q[0];
rz(-1.0210751) q[0];
sx q[0];
rz(1.8718029) q[0];
rz(-pi) q[1];
rz(-2.2631049) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(-2.0932587) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.6565894) q[1];
sx q[1];
rz(-2.2409229) q[1];
sx q[1];
rz(0.90799241) q[1];
x q[2];
rz(2.8273724) q[3];
sx q[3];
rz(-1.1952956) q[3];
sx q[3];
rz(-2.7995031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(-1.3683866) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(-0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1798582) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(2.7767048) q[0];
rz(-2.2015613) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.5023124) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2973328) q[0];
sx q[0];
rz(-1.5808006) q[0];
sx q[0];
rz(0.25015932) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83606007) q[2];
sx q[2];
rz(-1.2088472) q[2];
sx q[2];
rz(-1.8723633) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2545401) q[1];
sx q[1];
rz(-0.92003838) q[1];
sx q[1];
rz(-2.410143) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3650465) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(-2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6442948) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(2.0765182) q[2];
rz(-2.8403357) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97312462) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(1.3849974) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(2.7246144) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57396736) q[0];
sx q[0];
rz(-1.5331393) q[0];
sx q[0];
rz(2.2249939) q[0];
x q[1];
rz(0.68848227) q[2];
sx q[2];
rz(-0.79198972) q[2];
sx q[2];
rz(-0.85306963) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.13874395) q[1];
sx q[1];
rz(-1.4321623) q[1];
sx q[1];
rz(-2.0744051) q[1];
rz(3.0005089) q[3];
sx q[3];
rz(-1.699563) q[3];
sx q[3];
rz(-1.5878549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0372662) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(-2.9157675) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(-0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.726783) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(1.6171932) q[0];
rz(2.1879451) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(1.8189925) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1117301) q[0];
sx q[0];
rz(-2.7435281) q[0];
sx q[0];
rz(-1.8882621) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0595423) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(-0.79364712) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9433371) q[1];
sx q[1];
rz(-2.2554734) q[1];
sx q[1];
rz(0.68590045) q[1];
rz(2.2262276) q[3];
sx q[3];
rz(-1.4487106) q[3];
sx q[3];
rz(-1.9742427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(2.1255169) q[2];
rz(1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(-0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(2.6323742) q[2];
sx q[2];
rz(-1.5621395) q[2];
sx q[2];
rz(-0.12315673) q[2];
rz(1.8018467) q[3];
sx q[3];
rz(-1.6587202) q[3];
sx q[3];
rz(-0.29826577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
