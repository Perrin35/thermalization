OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0535799) q[0];
sx q[0];
rz(-0.3126643) q[0];
sx q[0];
rz(3.0076658) q[0];
rz(3.1349831) q[1];
sx q[1];
rz(-1.5899038) q[1];
sx q[1];
rz(-0.42491999) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6962587) q[0];
sx q[0];
rz(-0.15275341) q[0];
sx q[0];
rz(0.85102083) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0449465) q[2];
sx q[2];
rz(-2.4174848) q[2];
sx q[2];
rz(2.6158138) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1462789) q[1];
sx q[1];
rz(-1.9094719) q[1];
sx q[1];
rz(-0.92471735) q[1];
rz(-pi) q[2];
rz(2.4185804) q[3];
sx q[3];
rz(-2.0916846) q[3];
sx q[3];
rz(-0.82244825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6801844) q[2];
sx q[2];
rz(-0.37698656) q[2];
sx q[2];
rz(3.027463) q[2];
rz(2.7528609) q[3];
sx q[3];
rz(-0.26624334) q[3];
sx q[3];
rz(1.3346599) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1330133) q[0];
sx q[0];
rz(-2.7393434) q[0];
sx q[0];
rz(-0.026799686) q[0];
rz(1.9642448) q[1];
sx q[1];
rz(-0.64759308) q[1];
sx q[1];
rz(-0.037638232) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1397495) q[0];
sx q[0];
rz(-1.5842373) q[0];
sx q[0];
rz(0.77450181) q[0];
x q[1];
rz(-2.3956793) q[2];
sx q[2];
rz(-2.841171) q[2];
sx q[2];
rz(-2.3021509) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.85559713) q[1];
sx q[1];
rz(-1.0904795) q[1];
sx q[1];
rz(-2.2248286) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2977799) q[3];
sx q[3];
rz(-2.3721266) q[3];
sx q[3];
rz(0.96801341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8893975) q[2];
sx q[2];
rz(-2.1571721) q[2];
sx q[2];
rz(0.21749116) q[2];
rz(0.43976954) q[3];
sx q[3];
rz(-1.7871126) q[3];
sx q[3];
rz(-3.0396438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52461034) q[0];
sx q[0];
rz(-0.007402448) q[0];
sx q[0];
rz(-2.5819085) q[0];
rz(0.89645487) q[1];
sx q[1];
rz(-0.47210109) q[1];
sx q[1];
rz(2.7339973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85227634) q[0];
sx q[0];
rz(-1.6860776) q[0];
sx q[0];
rz(2.1570221) q[0];
rz(3.0733068) q[2];
sx q[2];
rz(-0.90203055) q[2];
sx q[2];
rz(-2.1055773) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7771908) q[1];
sx q[1];
rz(-0.47955105) q[1];
sx q[1];
rz(-0.24141356) q[1];
x q[2];
rz(-0.30826195) q[3];
sx q[3];
rz(-1.0438266) q[3];
sx q[3];
rz(-1.7402203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1233623) q[2];
sx q[2];
rz(-1.9841649) q[2];
sx q[2];
rz(2.1165712) q[2];
rz(1.7068663) q[3];
sx q[3];
rz(-0.44308174) q[3];
sx q[3];
rz(-0.6492492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0822815) q[0];
sx q[0];
rz(-2.870443) q[0];
sx q[0];
rz(-0.49736381) q[0];
rz(1.8479895) q[1];
sx q[1];
rz(-2.135364) q[1];
sx q[1];
rz(1.9088378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3022753) q[0];
sx q[0];
rz(-2.1324106) q[0];
sx q[0];
rz(1.8118493) q[0];
rz(-pi) q[1];
rz(-2.4380765) q[2];
sx q[2];
rz(-1.8986964) q[2];
sx q[2];
rz(-0.10645535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9962483) q[1];
sx q[1];
rz(-1.323902) q[1];
sx q[1];
rz(1.7801784) q[1];
rz(-pi) q[2];
rz(1.5895548) q[3];
sx q[3];
rz(-1.1705468) q[3];
sx q[3];
rz(1.7168604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5645912) q[2];
sx q[2];
rz(-0.61508721) q[2];
sx q[2];
rz(-2.9992529) q[2];
rz(1.9650991) q[3];
sx q[3];
rz(-2.1409972) q[3];
sx q[3];
rz(0.4308027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.38554) q[0];
sx q[0];
rz(-0.98243326) q[0];
sx q[0];
rz(2.8681712) q[0];
rz(-2.7418819) q[1];
sx q[1];
rz(-2.3276261) q[1];
sx q[1];
rz(0.25693691) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98824153) q[0];
sx q[0];
rz(-1.0437766) q[0];
sx q[0];
rz(-2.7668883) q[0];
rz(-0.24946282) q[2];
sx q[2];
rz(-1.1946861) q[2];
sx q[2];
rz(0.94406908) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.6338004) q[1];
sx q[1];
rz(-1.6163428) q[1];
sx q[1];
rz(-1.6409268) q[1];
rz(-pi) q[2];
rz(-2.2502348) q[3];
sx q[3];
rz(-1.7943503) q[3];
sx q[3];
rz(-0.67878658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5520681) q[2];
sx q[2];
rz(-0.35577154) q[2];
sx q[2];
rz(2.3815928) q[2];
rz(0.99203569) q[3];
sx q[3];
rz(-1.2090679) q[3];
sx q[3];
rz(-1.1442643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29042596) q[0];
sx q[0];
rz(-0.63153428) q[0];
sx q[0];
rz(2.5354711) q[0];
rz(1.0468613) q[1];
sx q[1];
rz(-1.4907587) q[1];
sx q[1];
rz(3.0643588) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.58409) q[0];
sx q[0];
rz(-1.5582005) q[0];
sx q[0];
rz(1.5303311) q[0];
rz(-pi) q[1];
rz(0.31862835) q[2];
sx q[2];
rz(-2.0214404) q[2];
sx q[2];
rz(-0.81540996) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97925767) q[1];
sx q[1];
rz(-1.1282776) q[1];
sx q[1];
rz(-2.6975432) q[1];
rz(-pi) q[2];
rz(2.4786754) q[3];
sx q[3];
rz(-1.7212711) q[3];
sx q[3];
rz(1.6126954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91654009) q[2];
sx q[2];
rz(-0.70987916) q[2];
sx q[2];
rz(0.27099657) q[2];
rz(0.58800507) q[3];
sx q[3];
rz(-2.3376412) q[3];
sx q[3];
rz(-0.40905455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8005017) q[0];
sx q[0];
rz(-1.611447) q[0];
sx q[0];
rz(-0.28277582) q[0];
rz(-0.68280363) q[1];
sx q[1];
rz(-2.2802201) q[1];
sx q[1];
rz(2.2883794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68180726) q[0];
sx q[0];
rz(-0.88621688) q[0];
sx q[0];
rz(2.6537978) q[0];
rz(-0.81099895) q[2];
sx q[2];
rz(-1.930522) q[2];
sx q[2];
rz(-0.0051509858) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68487271) q[1];
sx q[1];
rz(-1.3383075) q[1];
sx q[1];
rz(2.9732735) q[1];
rz(-pi) q[2];
rz(3.0158132) q[3];
sx q[3];
rz(-1.3759633) q[3];
sx q[3];
rz(1.3013713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0736488) q[2];
sx q[2];
rz(-2.6239519) q[2];
sx q[2];
rz(-0.29120564) q[2];
rz(1.9368885) q[3];
sx q[3];
rz(-2.5681345) q[3];
sx q[3];
rz(1.5706435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7912927) q[0];
sx q[0];
rz(-1.6181823) q[0];
sx q[0];
rz(0.6268025) q[0];
rz(2.0794012) q[1];
sx q[1];
rz(-0.79963446) q[1];
sx q[1];
rz(-0.83745426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27521321) q[0];
sx q[0];
rz(-1.4555706) q[0];
sx q[0];
rz(-2.000745) q[0];
rz(-pi) q[1];
rz(-0.59932389) q[2];
sx q[2];
rz(-2.7380652) q[2];
sx q[2];
rz(-2.7996792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9161842) q[1];
sx q[1];
rz(-1.7042158) q[1];
sx q[1];
rz(-2.0853985) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7193878) q[3];
sx q[3];
rz(-2.0839543) q[3];
sx q[3];
rz(-0.10529127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31967878) q[2];
sx q[2];
rz(-1.9210812) q[2];
sx q[2];
rz(2.1004045) q[2];
rz(-0.90211165) q[3];
sx q[3];
rz(-1.9769042) q[3];
sx q[3];
rz(0.56727099) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5015471) q[0];
sx q[0];
rz(-0.29842672) q[0];
sx q[0];
rz(-0.10064594) q[0];
rz(1.8177265) q[1];
sx q[1];
rz(-0.83815014) q[1];
sx q[1];
rz(0.59555882) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32257358) q[0];
sx q[0];
rz(-0.77054502) q[0];
sx q[0];
rz(0.84011232) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8331489) q[2];
sx q[2];
rz(-2.4624918) q[2];
sx q[2];
rz(-2.6237112) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2145077) q[1];
sx q[1];
rz(-0.86125285) q[1];
sx q[1];
rz(-1.2432712) q[1];
rz(-pi) q[2];
rz(1.4138537) q[3];
sx q[3];
rz(-2.3683386) q[3];
sx q[3];
rz(0.93676222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9969534) q[2];
sx q[2];
rz(-2.3914631) q[2];
sx q[2];
rz(-1.7993125) q[2];
rz(-0.73623776) q[3];
sx q[3];
rz(-2.8056371) q[3];
sx q[3];
rz(-3.0430702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1243684) q[0];
sx q[0];
rz(-2.4898744) q[0];
sx q[0];
rz(-2.8454054) q[0];
rz(-1.4172957) q[1];
sx q[1];
rz(-0.92767757) q[1];
sx q[1];
rz(-2.9790624) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2784341) q[0];
sx q[0];
rz(-1.5774283) q[0];
sx q[0];
rz(1.5445821) q[0];
x q[1];
rz(-2.1549781) q[2];
sx q[2];
rz(-0.59609813) q[2];
sx q[2];
rz(-1.7510027) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49419241) q[1];
sx q[1];
rz(-1.7583349) q[1];
sx q[1];
rz(1.4306551) q[1];
rz(-0.82270427) q[3];
sx q[3];
rz(-2.4260022) q[3];
sx q[3];
rz(-1.1038609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6235003) q[2];
sx q[2];
rz(-1.2178428) q[2];
sx q[2];
rz(-0.45563844) q[2];
rz(-1.0738922) q[3];
sx q[3];
rz(-0.62993252) q[3];
sx q[3];
rz(-0.58551252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.0475273) q[0];
sx q[0];
rz(-1.3205262) q[0];
sx q[0];
rz(2.5393215) q[0];
rz(2.8216254) q[1];
sx q[1];
rz(-2.0260369) q[1];
sx q[1];
rz(1.8021348) q[1];
rz(-1.6306277) q[2];
sx q[2];
rz(-1.9070121) q[2];
sx q[2];
rz(-0.4941834) q[2];
rz(1.6781758) q[3];
sx q[3];
rz(-0.79859514) q[3];
sx q[3];
rz(2.7664281) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
