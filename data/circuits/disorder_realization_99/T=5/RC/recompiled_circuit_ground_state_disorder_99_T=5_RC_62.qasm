OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91166624) q[0];
sx q[0];
rz(-0.32719964) q[0];
sx q[0];
rz(1.3258452) q[0];
rz(1.3658547) q[1];
sx q[1];
rz(-0.39443016) q[1];
sx q[1];
rz(-2.3529513) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3978826) q[0];
sx q[0];
rz(-1.2628444) q[0];
sx q[0];
rz(-0.90713199) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82420492) q[2];
sx q[2];
rz(-1.0366777) q[2];
sx q[2];
rz(-0.20301486) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4439075) q[1];
sx q[1];
rz(-0.34516066) q[1];
sx q[1];
rz(0.60072319) q[1];
x q[2];
rz(2.1988081) q[3];
sx q[3];
rz(-1.4664141) q[3];
sx q[3];
rz(-1.2918772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49393645) q[2];
sx q[2];
rz(-1.6948573) q[2];
sx q[2];
rz(-0.20230618) q[2];
rz(-0.10937396) q[3];
sx q[3];
rz(-2.5325363) q[3];
sx q[3];
rz(3.0561395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(-2.0789455) q[0];
sx q[0];
rz(-0.92342347) q[0];
sx q[0];
rz(1.1751291) q[0];
rz(-1.6773978) q[1];
sx q[1];
rz(-1.0025832) q[1];
sx q[1];
rz(2.7820803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4660366) q[0];
sx q[0];
rz(-0.17854843) q[0];
sx q[0];
rz(-1.1677062) q[0];
rz(1.1929252) q[2];
sx q[2];
rz(-2.0959575) q[2];
sx q[2];
rz(1.988387) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0766616) q[1];
sx q[1];
rz(-2.1691469) q[1];
sx q[1];
rz(-1.1538572) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1678108) q[3];
sx q[3];
rz(-1.0724073) q[3];
sx q[3];
rz(2.1417051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8466865) q[2];
sx q[2];
rz(-2.0945695) q[2];
sx q[2];
rz(-2.8666551) q[2];
rz(1.8168195) q[3];
sx q[3];
rz(-1.5271657) q[3];
sx q[3];
rz(1.2980609) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0637958) q[0];
sx q[0];
rz(-2.9834788) q[0];
sx q[0];
rz(-2.7445444) q[0];
rz(2.5322757) q[1];
sx q[1];
rz(-1.3343697) q[1];
sx q[1];
rz(1.5999925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55508321) q[0];
sx q[0];
rz(-1.6753046) q[0];
sx q[0];
rz(-1.2782607) q[0];
x q[1];
rz(-1.5055539) q[2];
sx q[2];
rz(-1.5582623) q[2];
sx q[2];
rz(1.6103075) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0761281) q[1];
sx q[1];
rz(-0.72678002) q[1];
sx q[1];
rz(1.5453592) q[1];
x q[2];
rz(-2.3261689) q[3];
sx q[3];
rz(-2.2924726) q[3];
sx q[3];
rz(1.4072925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93727532) q[2];
sx q[2];
rz(-1.5466362) q[2];
sx q[2];
rz(-0.23359648) q[2];
rz(1.8448081) q[3];
sx q[3];
rz(-1.1917043) q[3];
sx q[3];
rz(2.4209723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.63491708) q[0];
sx q[0];
rz(-1.9838061) q[0];
sx q[0];
rz(-1.3651715) q[0];
rz(-1.8695976) q[1];
sx q[1];
rz(-1.9422266) q[1];
sx q[1];
rz(-1.2619527) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1256169) q[0];
sx q[0];
rz(-1.7628094) q[0];
sx q[0];
rz(-0.58108347) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4883411) q[2];
sx q[2];
rz(-2.8318498) q[2];
sx q[2];
rz(2.7249634) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8820928) q[1];
sx q[1];
rz(-1.6354523) q[1];
sx q[1];
rz(-2.9586709) q[1];
rz(1.3586952) q[3];
sx q[3];
rz(-0.75772983) q[3];
sx q[3];
rz(1.084071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23123732) q[2];
sx q[2];
rz(-2.1114025) q[2];
sx q[2];
rz(2.3528698) q[2];
rz(1.8606868) q[3];
sx q[3];
rz(-1.3642045) q[3];
sx q[3];
rz(2.1196712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3874409) q[0];
sx q[0];
rz(-1.0197637) q[0];
sx q[0];
rz(-0.66106853) q[0];
rz(-2.0263653) q[1];
sx q[1];
rz(-1.3122357) q[1];
sx q[1];
rz(2.1818395) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5710053) q[0];
sx q[0];
rz(-2.5474605) q[0];
sx q[0];
rz(-2.6460365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8371253) q[2];
sx q[2];
rz(-2.5779471) q[2];
sx q[2];
rz(3.0802022) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.58127357) q[1];
sx q[1];
rz(-1.4380102) q[1];
sx q[1];
rz(1.8627235) q[1];
x q[2];
rz(1.3101391) q[3];
sx q[3];
rz(-0.66124004) q[3];
sx q[3];
rz(-0.52296296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5614718) q[2];
sx q[2];
rz(-0.46130195) q[2];
sx q[2];
rz(2.0816154) q[2];
rz(-2.3914242) q[3];
sx q[3];
rz(-1.7629905) q[3];
sx q[3];
rz(1.2257956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8691413) q[0];
sx q[0];
rz(-1.5584109) q[0];
sx q[0];
rz(2.9172752) q[0];
rz(-1.6250826) q[1];
sx q[1];
rz(-2.5254011) q[1];
sx q[1];
rz(3.065965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6859378) q[0];
sx q[0];
rz(-1.9561531) q[0];
sx q[0];
rz(-2.1598201) q[0];
rz(-pi) q[1];
rz(-1.8746156) q[2];
sx q[2];
rz(-1.4785789) q[2];
sx q[2];
rz(0.95668787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4757753) q[1];
sx q[1];
rz(-0.66002895) q[1];
sx q[1];
rz(2.3296861) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7073955) q[3];
sx q[3];
rz(-0.89354529) q[3];
sx q[3];
rz(1.7627258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.533796) q[2];
sx q[2];
rz(-0.7527315) q[2];
sx q[2];
rz(0.047018615) q[2];
rz(2.4646344) q[3];
sx q[3];
rz(-1.2983863) q[3];
sx q[3];
rz(0.025378749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(0.15636477) q[0];
sx q[0];
rz(-1.6522836) q[0];
sx q[0];
rz(1.697668) q[0];
rz(-2.2185183) q[1];
sx q[1];
rz(-2.1363027) q[1];
sx q[1];
rz(2.9583171) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61675894) q[0];
sx q[0];
rz(-1.3206894) q[0];
sx q[0];
rz(2.564774) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44281339) q[2];
sx q[2];
rz(-1.9573136) q[2];
sx q[2];
rz(2.7903976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.029327) q[1];
sx q[1];
rz(-1.8958175) q[1];
sx q[1];
rz(2.0954262) q[1];
rz(-2.8541862) q[3];
sx q[3];
rz(-1.0106405) q[3];
sx q[3];
rz(2.0623824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1075333) q[2];
sx q[2];
rz(-1.6851765) q[2];
sx q[2];
rz(0.019406645) q[2];
rz(0.16128811) q[3];
sx q[3];
rz(-1.015377) q[3];
sx q[3];
rz(1.2279145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0374544) q[0];
sx q[0];
rz(-0.4087953) q[0];
sx q[0];
rz(3.0492875) q[0];
rz(-1.1760938) q[1];
sx q[1];
rz(-0.25830019) q[1];
sx q[1];
rz(1.7324804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8831439) q[0];
sx q[0];
rz(-1.5165268) q[0];
sx q[0];
rz(-1.4656214) q[0];
rz(0.46940501) q[2];
sx q[2];
rz(-1.0115716) q[2];
sx q[2];
rz(2.6162868) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8586052) q[1];
sx q[1];
rz(-1.8388565) q[1];
sx q[1];
rz(0.40045935) q[1];
rz(-pi) q[2];
rz(-2.4020089) q[3];
sx q[3];
rz(-0.7639262) q[3];
sx q[3];
rz(2.6134864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.041542355) q[2];
sx q[2];
rz(-1.5644194) q[2];
sx q[2];
rz(-1.9680295) q[2];
rz(0.38343492) q[3];
sx q[3];
rz(-1.8337367) q[3];
sx q[3];
rz(-2.0508544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1293056) q[0];
sx q[0];
rz(-1.4947991) q[0];
sx q[0];
rz(2.9360085) q[0];
rz(2.3221817) q[1];
sx q[1];
rz(-1.2148379) q[1];
sx q[1];
rz(-2.6507071) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48614472) q[0];
sx q[0];
rz(-1.9441588) q[0];
sx q[0];
rz(0.70151897) q[0];
rz(-pi) q[1];
rz(1.8794291) q[2];
sx q[2];
rz(-2.1763945) q[2];
sx q[2];
rz(-0.016589368) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1236022) q[1];
sx q[1];
rz(-1.523087) q[1];
sx q[1];
rz(-2.1393659) q[1];
x q[2];
rz(1.396046) q[3];
sx q[3];
rz(-0.51769231) q[3];
sx q[3];
rz(-0.92624901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72863355) q[2];
sx q[2];
rz(-2.1042991) q[2];
sx q[2];
rz(1.884985) q[2];
rz(2.7058153) q[3];
sx q[3];
rz(-1.1083009) q[3];
sx q[3];
rz(-0.88424879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7213781) q[0];
sx q[0];
rz(-0.90737897) q[0];
sx q[0];
rz(-0.23289982) q[0];
rz(0.13889343) q[1];
sx q[1];
rz(-1.9878309) q[1];
sx q[1];
rz(-0.5084261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4111944) q[0];
sx q[0];
rz(-1.8286213) q[0];
sx q[0];
rz(1.0143552) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1022212) q[2];
sx q[2];
rz(-0.55742427) q[2];
sx q[2];
rz(-2.1845326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99791807) q[1];
sx q[1];
rz(-1.7047722) q[1];
sx q[1];
rz(-0.0079946144) q[1];
rz(-0.53896972) q[3];
sx q[3];
rz(-2.1930755) q[3];
sx q[3];
rz(1.8405746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.4109219) q[2];
sx q[2];
rz(-1.4620917) q[2];
sx q[2];
rz(-0.81653583) q[2];
rz(0.38980347) q[3];
sx q[3];
rz(-2.1392418) q[3];
sx q[3];
rz(-1.3780814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941866) q[0];
sx q[0];
rz(-2.2185855) q[0];
sx q[0];
rz(2.9119281) q[0];
rz(-1.6387088) q[1];
sx q[1];
rz(-1.8062183) q[1];
sx q[1];
rz(0.91581215) q[1];
rz(2.7336929) q[2];
sx q[2];
rz(-0.41389155) q[2];
sx q[2];
rz(2.2663739) q[2];
rz(1.3744204) q[3];
sx q[3];
rz(-2.4300601) q[3];
sx q[3];
rz(-2.0236494) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
