OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52743071) q[0];
sx q[0];
rz(3.951374) q[0];
sx q[0];
rz(9.9561719) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(-1.7778492) q[1];
sx q[1];
rz(-1.2030503) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1019842) q[0];
sx q[0];
rz(-1.7829722) q[0];
sx q[0];
rz(-0.19594812) q[0];
x q[1];
rz(-1.7877903) q[2];
sx q[2];
rz(-1.8197682) q[2];
sx q[2];
rz(2.4169902) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4431475) q[1];
sx q[1];
rz(-1.5138211) q[1];
sx q[1];
rz(0.2351825) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0193127) q[3];
sx q[3];
rz(-0.29502007) q[3];
sx q[3];
rz(2.2515841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8649341) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(2.9425088) q[2];
rz(-1.52786) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(0.73408192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1428225) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(0.81940991) q[0];
rz(-2.8557414) q[1];
sx q[1];
rz(-0.91393036) q[1];
sx q[1];
rz(1.8751289) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0430543) q[0];
sx q[0];
rz(-2.1715954) q[0];
sx q[0];
rz(-2.0354802) q[0];
rz(1.0533633) q[2];
sx q[2];
rz(-0.45767637) q[2];
sx q[2];
rz(-0.44331726) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4283735) q[1];
sx q[1];
rz(-1.1717147) q[1];
sx q[1];
rz(0.39308163) q[1];
x q[2];
rz(0.5760848) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(1.9379804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1566029) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(1.161969) q[2];
rz(-3.0544288) q[3];
sx q[3];
rz(-1.6379084) q[3];
sx q[3];
rz(-0.073908977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6067628) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(0.30360046) q[0];
rz(-1.7595694) q[1];
sx q[1];
rz(-1.8654114) q[1];
sx q[1];
rz(0.64750013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.835639) q[0];
sx q[0];
rz(-1.5257972) q[0];
sx q[0];
rz(0.63458981) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60909231) q[2];
sx q[2];
rz(-0.81431544) q[2];
sx q[2];
rz(2.7012205) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4422569) q[1];
sx q[1];
rz(-1.934821) q[1];
sx q[1];
rz(2.060021) q[1];
rz(-pi) q[2];
rz(-2.2494499) q[3];
sx q[3];
rz(-2.3193079) q[3];
sx q[3];
rz(-0.014978623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8422164) q[2];
sx q[2];
rz(-0.78263485) q[2];
sx q[2];
rz(1.5608609) q[2];
rz(-2.1598699) q[3];
sx q[3];
rz(-1.862062) q[3];
sx q[3];
rz(-2.4981807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2340387) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(-2.2696944) q[0];
rz(-2.7543228) q[1];
sx q[1];
rz(-0.64278066) q[1];
sx q[1];
rz(0.29104582) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5975808) q[0];
sx q[0];
rz(-0.62950069) q[0];
sx q[0];
rz(-1.9660216) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8918858) q[2];
sx q[2];
rz(-1.849527) q[2];
sx q[2];
rz(-2.2433787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9705022) q[1];
sx q[1];
rz(-0.88741747) q[1];
sx q[1];
rz(0.6716397) q[1];
rz(-1.325874) q[3];
sx q[3];
rz(-2.5513253) q[3];
sx q[3];
rz(3.0274689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7724093) q[2];
sx q[2];
rz(-0.79493752) q[2];
sx q[2];
rz(-2.5975361) q[2];
rz(-2.6323075) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(-2.2756186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13957025) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(0.64506662) q[0];
rz(0.6257261) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(-2.5114139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1030578) q[0];
sx q[0];
rz(-1.6928506) q[0];
sx q[0];
rz(2.5621668) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9899881) q[2];
sx q[2];
rz(-2.345511) q[2];
sx q[2];
rz(2.5615356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.899257) q[1];
sx q[1];
rz(-1.7944272) q[1];
sx q[1];
rz(1.0244589) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24124055) q[3];
sx q[3];
rz(-2.5706228) q[3];
sx q[3];
rz(0.2824479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75382918) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(-1.3641664) q[2];
rz(-1.8917313) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(0.92145872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0867778) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(-0.71682799) q[0];
rz(2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(2.3278918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98052927) q[0];
sx q[0];
rz(-1.5080394) q[0];
sx q[0];
rz(-1.784523) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6394303) q[2];
sx q[2];
rz(-0.8317906) q[2];
sx q[2];
rz(-0.30314988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54436436) q[1];
sx q[1];
rz(-1.2372984) q[1];
sx q[1];
rz(0.79798214) q[1];
rz(-pi) q[2];
rz(-2.8343809) q[3];
sx q[3];
rz(-1.2324411) q[3];
sx q[3];
rz(-0.070778155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2563236) q[2];
sx q[2];
rz(-0.43748125) q[2];
sx q[2];
rz(-1.9539333) q[2];
rz(0.93196431) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1535783) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(2.4555092) q[0];
rz(-1.5230806) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(0.66326052) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30090573) q[0];
sx q[0];
rz(-1.5468883) q[0];
sx q[0];
rz(-0.92047128) q[0];
rz(-pi) q[1];
rz(2.5040607) q[2];
sx q[2];
rz(-1.3946643) q[2];
sx q[2];
rz(-2.9498364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82603964) q[1];
sx q[1];
rz(-1.3693046) q[1];
sx q[1];
rz(1.0793346) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0955986) q[3];
sx q[3];
rz(-1.0605269) q[3];
sx q[3];
rz(1.0218395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66701165) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(-2.1298501) q[3];
sx q[3];
rz(-1.116131) q[3];
sx q[3];
rz(-2.570178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0886154) q[0];
sx q[0];
rz(-2.4089854) q[0];
sx q[0];
rz(0.10738871) q[0];
rz(2.836851) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(-2.6838141) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8900507) q[0];
sx q[0];
rz(-0.72351447) q[0];
sx q[0];
rz(1.0735805) q[0];
rz(-2.3950251) q[2];
sx q[2];
rz(-1.1636359) q[2];
sx q[2];
rz(-2.6395869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0486054) q[1];
sx q[1];
rz(-2.7106206) q[1];
sx q[1];
rz(2.155817) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9321953) q[3];
sx q[3];
rz(-2.0933588) q[3];
sx q[3];
rz(0.88025974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7765939) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.4314852) q[2];
rz(1.7163904) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9872221) q[0];
sx q[0];
rz(-0.4796589) q[0];
sx q[0];
rz(2.1881058) q[0];
rz(1.8158688) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(0.58473933) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3922694) q[0];
sx q[0];
rz(-1.4587147) q[0];
sx q[0];
rz(-2.6507676) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5647613) q[2];
sx q[2];
rz(-0.62099651) q[2];
sx q[2];
rz(-0.15658212) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.016420267) q[1];
sx q[1];
rz(-2.6965045) q[1];
sx q[1];
rz(-0.16146407) q[1];
rz(-pi) q[2];
rz(2.0546012) q[3];
sx q[3];
rz(-0.30011794) q[3];
sx q[3];
rz(0.3008315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3328302) q[2];
sx q[2];
rz(-2.3213883) q[2];
sx q[2];
rz(2.8653223) q[2];
rz(2.590498) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(2.7579257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0693414) q[0];
sx q[0];
rz(-2.6273917) q[0];
sx q[0];
rz(2.4861091) q[0];
rz(1.1570702) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(1.5225333) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039278395) q[0];
sx q[0];
rz(-0.73844665) q[0];
sx q[0];
rz(0.46955903) q[0];
rz(-1.9117781) q[2];
sx q[2];
rz(-0.91170646) q[2];
sx q[2];
rz(1.9552719) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9177502) q[1];
sx q[1];
rz(-2.192944) q[1];
sx q[1];
rz(1.7734581) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7195815) q[3];
sx q[3];
rz(-1.2000298) q[3];
sx q[3];
rz(-1.5164204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15554252) q[2];
sx q[2];
rz(-0.4077929) q[2];
sx q[2];
rz(-0.84214169) q[2];
rz(-1.5367674) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(3.1150637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90606541) q[0];
sx q[0];
rz(-1.9530095) q[0];
sx q[0];
rz(-0.45146913) q[0];
rz(0.72262598) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(1.6956383) q[2];
sx q[2];
rz(-2.9308133) q[2];
sx q[2];
rz(1.4598893) q[2];
rz(1.6668672) q[3];
sx q[3];
rz(-2.4484652) q[3];
sx q[3];
rz(-0.0948003) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
