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
rz(1.2713852) q[0];
sx q[0];
rz(-0.013590824) q[0];
sx q[0];
rz(3.115227) q[0];
rz(0.68323505) q[1];
sx q[1];
rz(-1.3807715) q[1];
sx q[1];
rz(0.071579054) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1522823) q[0];
sx q[0];
rz(-1.5805065) q[0];
sx q[0];
rz(-3.1165313) q[0];
rz(-pi) q[1];
rz(1.8060922) q[2];
sx q[2];
rz(-1.9517731) q[2];
sx q[2];
rz(1.345696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84715473) q[1];
sx q[1];
rz(-1.5757474) q[1];
sx q[1];
rz(1.5895542) q[1];
rz(-1.2504225) q[3];
sx q[3];
rz(-1.2696869) q[3];
sx q[3];
rz(1.4442577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2804395) q[2];
sx q[2];
rz(-0.70715487) q[2];
sx q[2];
rz(-0.46671483) q[2];
rz(-0.47433445) q[3];
sx q[3];
rz(-0.021183906) q[3];
sx q[3];
rz(3.1202313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83653432) q[0];
sx q[0];
rz(-2.6439522) q[0];
sx q[0];
rz(-3.1217788) q[0];
rz(-1.5927947) q[1];
sx q[1];
rz(-2.9217547) q[1];
sx q[1];
rz(1.4733431) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1602185) q[0];
sx q[0];
rz(-2.2607973) q[0];
sx q[0];
rz(-0.86440683) q[0];
x q[1];
rz(1.8657487) q[2];
sx q[2];
rz(-1.3222857) q[2];
sx q[2];
rz(-2.4139443) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8791318) q[1];
sx q[1];
rz(-1.3768798) q[1];
sx q[1];
rz(0.078207774) q[1];
x q[2];
rz(-0.83062828) q[3];
sx q[3];
rz(-2.1601956) q[3];
sx q[3];
rz(-0.43260655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8389429) q[2];
sx q[2];
rz(-0.5984211) q[2];
sx q[2];
rz(-1.8518651) q[2];
rz(1.2368115) q[3];
sx q[3];
rz(-0.33331063) q[3];
sx q[3];
rz(2.439177) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1692093) q[0];
sx q[0];
rz(-1.9820259) q[0];
sx q[0];
rz(1.5709391) q[0];
rz(1.6698569) q[1];
sx q[1];
rz(-1.6153299) q[1];
sx q[1];
rz(2.6872046) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4226625) q[0];
sx q[0];
rz(-0.61394982) q[0];
sx q[0];
rz(-1.802868) q[0];
rz(-pi) q[1];
x q[1];
rz(1.387792) q[2];
sx q[2];
rz(-2.9972074) q[2];
sx q[2];
rz(-0.78156495) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4207698) q[1];
sx q[1];
rz(-1.5934363) q[1];
sx q[1];
rz(-1.6890425) q[1];
rz(0.91505312) q[3];
sx q[3];
rz(-2.1106535) q[3];
sx q[3];
rz(1.8600382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4001974) q[2];
sx q[2];
rz(-0.016308451) q[2];
sx q[2];
rz(2.7941217) q[2];
rz(0.45626429) q[3];
sx q[3];
rz(-0.01472344) q[3];
sx q[3];
rz(-0.99304503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72661138) q[0];
sx q[0];
rz(-1.9009637) q[0];
sx q[0];
rz(-1.4062784) q[0];
rz(2.6982488) q[1];
sx q[1];
rz(-2.1163546) q[1];
sx q[1];
rz(1.5700856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8271874) q[0];
sx q[0];
rz(-1.5159722) q[0];
sx q[0];
rz(-2.501295) q[0];
rz(-pi) q[1];
rz(-2.3894044) q[2];
sx q[2];
rz(-3.0471932) q[2];
sx q[2];
rz(-2.029325) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0519331) q[1];
sx q[1];
rz(-1.5610236) q[1];
sx q[1];
rz(-1.8490514) q[1];
rz(-pi) q[2];
rz(2.8624503) q[3];
sx q[3];
rz(-1.1542392) q[3];
sx q[3];
rz(-0.032241658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7683679) q[2];
sx q[2];
rz(-2.7478605) q[2];
sx q[2];
rz(0.079252871) q[2];
rz(-1.9521889) q[3];
sx q[3];
rz(-1.3711843) q[3];
sx q[3];
rz(-1.6325379) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4382512) q[0];
sx q[0];
rz(-0.51333135) q[0];
sx q[0];
rz(2.3355423) q[0];
rz(0.84000677) q[1];
sx q[1];
rz(-0.012931074) q[1];
sx q[1];
rz(-0.77617019) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66444976) q[0];
sx q[0];
rz(-2.5390671) q[0];
sx q[0];
rz(2.9092034) q[0];
rz(0.1397392) q[2];
sx q[2];
rz(-3.1310509) q[2];
sx q[2];
rz(0.14103061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8373903) q[1];
sx q[1];
rz(-1.5826962) q[1];
sx q[1];
rz(3.0061327) q[1];
rz(-pi) q[2];
rz(0.25960323) q[3];
sx q[3];
rz(-1.6791376) q[3];
sx q[3];
rz(-1.3110127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.025803056) q[2];
sx q[2];
rz(-1.5805406) q[2];
sx q[2];
rz(2.4157794) q[2];
rz(-0.18802655) q[3];
sx q[3];
rz(-3.0811716) q[3];
sx q[3];
rz(-2.4316725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96227729) q[0];
sx q[0];
rz(-0.55914068) q[0];
sx q[0];
rz(-0.54139262) q[0];
rz(2.9560282) q[1];
sx q[1];
rz(-1.5488397) q[1];
sx q[1];
rz(-3.013179) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5207883) q[0];
sx q[0];
rz(-2.6628011) q[0];
sx q[0];
rz(2.4420681) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4494684) q[2];
sx q[2];
rz(-1.5660966) q[2];
sx q[2];
rz(0.0041065816) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7306108) q[1];
sx q[1];
rz(-1.990098) q[1];
sx q[1];
rz(-1.6602519) q[1];
rz(-1.6261934) q[3];
sx q[3];
rz(-2.9528816) q[3];
sx q[3];
rz(2.1907937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7554756) q[2];
sx q[2];
rz(-3.0839034) q[2];
sx q[2];
rz(2.3003787) q[2];
rz(-0.20127131) q[3];
sx q[3];
rz(-1.6095251) q[3];
sx q[3];
rz(2.8819528) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5614618) q[0];
sx q[0];
rz(-0.78829563) q[0];
sx q[0];
rz(1.5808251) q[0];
rz(0.67281094) q[1];
sx q[1];
rz(-1.7377868) q[1];
sx q[1];
rz(0.028884551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0085379) q[0];
sx q[0];
rz(-1.2875746) q[0];
sx q[0];
rz(1.828156) q[0];
rz(1.3703652) q[2];
sx q[2];
rz(-1.1005963) q[2];
sx q[2];
rz(-1.4154132) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1950732) q[1];
sx q[1];
rz(-1.9819489) q[1];
sx q[1];
rz(2.8821936) q[1];
rz(-pi) q[2];
rz(-0.77856346) q[3];
sx q[3];
rz(-1.9664008) q[3];
sx q[3];
rz(-2.9958519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9081356) q[2];
sx q[2];
rz(-0.59671777) q[2];
sx q[2];
rz(2.2133568) q[2];
rz(2.8440031) q[3];
sx q[3];
rz(-2.987515) q[3];
sx q[3];
rz(0.84460622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069227844) q[0];
sx q[0];
rz(-2.897825) q[0];
sx q[0];
rz(-0.099076554) q[0];
rz(-2.15436) q[1];
sx q[1];
rz(-1.3039373) q[1];
sx q[1];
rz(2.6027021) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9469556) q[0];
sx q[0];
rz(-1.6536923) q[0];
sx q[0];
rz(-0.006587365) q[0];
rz(-pi) q[1];
rz(-2.0343696) q[2];
sx q[2];
rz(-0.055214334) q[2];
sx q[2];
rz(-2.5718073) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.423279) q[1];
sx q[1];
rz(-2.4850902) q[1];
sx q[1];
rz(2.5274407) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5032926) q[3];
sx q[3];
rz(-1.5482404) q[3];
sx q[3];
rz(1.8543138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.146356) q[2];
sx q[2];
rz(-3.1260999) q[2];
sx q[2];
rz(2.8323979) q[2];
rz(-2.501798) q[3];
sx q[3];
rz(-0.00034172405) q[3];
sx q[3];
rz(0.063808002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30534202) q[0];
sx q[0];
rz(-0.59665614) q[0];
sx q[0];
rz(0.054656595) q[0];
rz(2.0089741) q[1];
sx q[1];
rz(-1.9839958) q[1];
sx q[1];
rz(1.8554912) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4187546) q[0];
sx q[0];
rz(-1.6143454) q[0];
sx q[0];
rz(1.7489793) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80100061) q[2];
sx q[2];
rz(-3.082976) q[2];
sx q[2];
rz(0.80551565) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9735582) q[1];
sx q[1];
rz(-1.66093) q[1];
sx q[1];
rz(2.9154577) q[1];
rz(-pi) q[2];
rz(2.0580106) q[3];
sx q[3];
rz(-1.5069252) q[3];
sx q[3];
rz(-2.4364249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9201811) q[2];
sx q[2];
rz(-2.5765918) q[2];
sx q[2];
rz(-2.0378225) q[2];
rz(-1.53299) q[3];
sx q[3];
rz(-0.04403232) q[3];
sx q[3];
rz(-2.485763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-3.1405545) q[0];
sx q[0];
rz(-2.9619205) q[0];
sx q[0];
rz(3.1371064) q[0];
rz(-1.5525612) q[1];
sx q[1];
rz(-1.6922502) q[1];
sx q[1];
rz(-3.0844614) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4085081) q[0];
sx q[0];
rz(-1.391469) q[0];
sx q[0];
rz(-0.097436949) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9190019) q[2];
sx q[2];
rz(-1.711297) q[2];
sx q[2];
rz(-0.3046356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2559214) q[1];
sx q[1];
rz(-0.90270611) q[1];
sx q[1];
rz(0.1541962) q[1];
rz(-pi) q[2];
rz(1.7369676) q[3];
sx q[3];
rz(-2.3944085) q[3];
sx q[3];
rz(0.93624828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39157465) q[2];
sx q[2];
rz(-0.027316814) q[2];
sx q[2];
rz(-2.2726783) q[2];
rz(2.1598375) q[3];
sx q[3];
rz(-0.029564094) q[3];
sx q[3];
rz(2.6386007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045573087) q[0];
sx q[0];
rz(-1.4843142) q[0];
sx q[0];
rz(1.6577161) q[0];
rz(-0.45687301) q[1];
sx q[1];
rz(-0.15468205) q[1];
sx q[1];
rz(-0.044943132) q[1];
rz(-2.9484684) q[2];
sx q[2];
rz(-0.77342351) q[2];
sx q[2];
rz(0.20062994) q[2];
rz(1.7140688) q[3];
sx q[3];
rz(-2.2289056) q[3];
sx q[3];
rz(-3.0994305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
