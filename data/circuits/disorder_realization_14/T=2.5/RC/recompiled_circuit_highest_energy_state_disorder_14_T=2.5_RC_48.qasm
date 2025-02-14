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
rz(2.4013588) q[0];
sx q[0];
rz(-1.6594247) q[0];
sx q[0];
rz(0.33492127) q[0];
rz(0.51796335) q[1];
sx q[1];
rz(5.2809102) q[1];
sx q[1];
rz(10.032293) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018463919) q[0];
sx q[0];
rz(-1.0842609) q[0];
sx q[0];
rz(2.8790265) q[0];
x q[1];
rz(1.363475) q[2];
sx q[2];
rz(-0.51864114) q[2];
sx q[2];
rz(-0.27388369) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.78215296) q[1];
sx q[1];
rz(-2.6342168) q[1];
sx q[1];
rz(-0.33276593) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0161765) q[3];
sx q[3];
rz(-1.2313255) q[3];
sx q[3];
rz(2.826863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6174378) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(-0.66317916) q[2];
rz(-3.0607306) q[3];
sx q[3];
rz(-0.20564779) q[3];
sx q[3];
rz(-1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8993768) q[0];
sx q[0];
rz(-1.3726534) q[0];
sx q[0];
rz(-0.85897613) q[0];
rz(1.2902749) q[1];
sx q[1];
rz(-1.4651508) q[1];
sx q[1];
rz(-1.4069517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5840184) q[0];
sx q[0];
rz(-1.0913335) q[0];
sx q[0];
rz(0.25734253) q[0];
rz(0.036769899) q[2];
sx q[2];
rz(-1.2388133) q[2];
sx q[2];
rz(2.9197502) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21995658) q[1];
sx q[1];
rz(-2.0825279) q[1];
sx q[1];
rz(-0.10351609) q[1];
rz(-pi) q[2];
rz(-0.59515335) q[3];
sx q[3];
rz(-2.4754324) q[3];
sx q[3];
rz(2.7056138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0591639) q[2];
sx q[2];
rz(-2.6291206) q[2];
sx q[2];
rz(-1.814369) q[2];
rz(2.0969157) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(0.41675848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.5752207) q[0];
sx q[0];
rz(-0.91082585) q[0];
sx q[0];
rz(2.7161993) q[0];
rz(-1.7644024) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(-1.4345217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1578428) q[0];
sx q[0];
rz(-2.2089777) q[0];
sx q[0];
rz(-1.7179836) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6041669) q[2];
sx q[2];
rz(-0.86276189) q[2];
sx q[2];
rz(1.58085) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.94043865) q[1];
sx q[1];
rz(-2.33719) q[1];
sx q[1];
rz(-2.5522347) q[1];
x q[2];
rz(1.9886964) q[3];
sx q[3];
rz(-2.4800081) q[3];
sx q[3];
rz(-2.083287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44125685) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(-1.2909935) q[2];
rz(1.357632) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62035471) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(-0.77019101) q[0];
rz(0.99984804) q[1];
sx q[1];
rz(-2.5412173) q[1];
sx q[1];
rz(-1.8005449) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0324381) q[0];
sx q[0];
rz(-2.1103835) q[0];
sx q[0];
rz(1.3225609) q[0];
x q[1];
rz(-0.45939873) q[2];
sx q[2];
rz(-2.2041577) q[2];
sx q[2];
rz(-0.49672302) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8902258) q[1];
sx q[1];
rz(-2.736675) q[1];
sx q[1];
rz(0.7271073) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2854659) q[3];
sx q[3];
rz(-1.6661246) q[3];
sx q[3];
rz(-2.8872629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8009214) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(0.7684024) q[2];
rz(-3.0411804) q[3];
sx q[3];
rz(-1.6008335) q[3];
sx q[3];
rz(1.8628619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95555821) q[0];
sx q[0];
rz(-0.30534196) q[0];
sx q[0];
rz(0.22698639) q[0];
rz(1.3817894) q[1];
sx q[1];
rz(-2.558936) q[1];
sx q[1];
rz(-1.3963799) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5433301) q[0];
sx q[0];
rz(-0.51307438) q[0];
sx q[0];
rz(-2.3725879) q[0];
rz(-pi) q[1];
rz(-1.7368083) q[2];
sx q[2];
rz(-1.1653656) q[2];
sx q[2];
rz(2.977598) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.13670838) q[1];
sx q[1];
rz(-1.6423823) q[1];
sx q[1];
rz(0.63386713) q[1];
x q[2];
rz(-1.5377858) q[3];
sx q[3];
rz(-2.3950658) q[3];
sx q[3];
rz(-1.3804264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9153626) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(-3.0653595) q[2];
rz(-1.3876312) q[3];
sx q[3];
rz(-0.97532719) q[3];
sx q[3];
rz(-1.7414198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48080322) q[0];
sx q[0];
rz(-1.5605518) q[0];
sx q[0];
rz(-1.3355108) q[0];
rz(-0.90323365) q[1];
sx q[1];
rz(-1.8025554) q[1];
sx q[1];
rz(-0.18641557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6680697) q[0];
sx q[0];
rz(-2.0041564) q[0];
sx q[0];
rz(-2.1587203) q[0];
rz(0.30419402) q[2];
sx q[2];
rz(-2.549571) q[2];
sx q[2];
rz(-1.5409868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32115667) q[1];
sx q[1];
rz(-1.8119436) q[1];
sx q[1];
rz(-2.2187869) q[1];
x q[2];
rz(2.0669492) q[3];
sx q[3];
rz(-2.0130139) q[3];
sx q[3];
rz(-2.7520455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4621801) q[2];
sx q[2];
rz(-1.1932411) q[2];
sx q[2];
rz(-1.818044) q[2];
rz(1.9994252) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(-2.2511258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46397504) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(0.63419813) q[0];
rz(2.0967261) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(-0.96010906) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30859892) q[0];
sx q[0];
rz(-1.9694298) q[0];
sx q[0];
rz(0.029775814) q[0];
rz(-pi) q[1];
rz(0.30107408) q[2];
sx q[2];
rz(-1.5784401) q[2];
sx q[2];
rz(-2.6507225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3297616) q[1];
sx q[1];
rz(-1.0454181) q[1];
sx q[1];
rz(3.1282022) q[1];
rz(-pi) q[2];
rz(-0.76310254) q[3];
sx q[3];
rz(-2.1402485) q[3];
sx q[3];
rz(0.72009898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.94379696) q[2];
sx q[2];
rz(-2.4891977) q[2];
sx q[2];
rz(1.5956399) q[2];
rz(-1.724285) q[3];
sx q[3];
rz(-0.82481074) q[3];
sx q[3];
rz(1.6212757) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6087795) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(2.1667495) q[0];
rz(0.10803647) q[1];
sx q[1];
rz(-1.255475) q[1];
sx q[1];
rz(1.1688165) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2377396) q[0];
sx q[0];
rz(-1.0971591) q[0];
sx q[0];
rz(2.6698378) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.37974) q[2];
sx q[2];
rz(-1.5708062) q[2];
sx q[2];
rz(2.5229483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.16202422) q[1];
sx q[1];
rz(-1.8037533) q[1];
sx q[1];
rz(-3.0515025) q[1];
rz(1.5904558) q[3];
sx q[3];
rz(-1.0791313) q[3];
sx q[3];
rz(1.8976854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.03269) q[2];
sx q[2];
rz(-2.2638075) q[2];
sx q[2];
rz(0.6558134) q[2];
rz(-3.0374895) q[3];
sx q[3];
rz(-1.7963573) q[3];
sx q[3];
rz(-0.18812215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26466894) q[0];
sx q[0];
rz(-1.8380565) q[0];
sx q[0];
rz(-1.4917829) q[0];
rz(-1.0393633) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(2.6731491) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3062564) q[0];
sx q[0];
rz(-0.1876615) q[0];
sx q[0];
rz(3.1188278) q[0];
x q[1];
rz(2.5889187) q[2];
sx q[2];
rz(-0.50853339) q[2];
sx q[2];
rz(2.7486211) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0026928) q[1];
sx q[1];
rz(-1.4093834) q[1];
sx q[1];
rz(2.3494865) q[1];
x q[2];
rz(-1.8519206) q[3];
sx q[3];
rz(-1.8418962) q[3];
sx q[3];
rz(-0.93096126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.09482065) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(-0.92791933) q[2];
rz(1.3377442) q[3];
sx q[3];
rz(-0.51023054) q[3];
sx q[3];
rz(1.4891589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054166404) q[0];
sx q[0];
rz(-1.0324284) q[0];
sx q[0];
rz(-2.0637276) q[0];
rz(2.7742591) q[1];
sx q[1];
rz(-1.8108188) q[1];
sx q[1];
rz(1.2841388) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3340942) q[0];
sx q[0];
rz(-2.0413412) q[0];
sx q[0];
rz(-2.8583291) q[0];
rz(1.7214891) q[2];
sx q[2];
rz(-0.65905276) q[2];
sx q[2];
rz(-2.6113752) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2312647) q[1];
sx q[1];
rz(-0.18016768) q[1];
sx q[1];
rz(2.641201) q[1];
x q[2];
rz(-1.7824836) q[3];
sx q[3];
rz(-0.84722391) q[3];
sx q[3];
rz(2.7016957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1401356) q[2];
sx q[2];
rz(-0.45150253) q[2];
sx q[2];
rz(2.7122811) q[2];
rz(-2.9863827) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(1.3252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975288) q[0];
sx q[0];
rz(-1.5915992) q[0];
sx q[0];
rz(-1.5912548) q[0];
rz(2.2776729) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(1.3270039) q[2];
sx q[2];
rz(-2.1012123) q[2];
sx q[2];
rz(-0.13441555) q[2];
rz(-2.8020017) q[3];
sx q[3];
rz(-2.17008) q[3];
sx q[3];
rz(-1.9934987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
