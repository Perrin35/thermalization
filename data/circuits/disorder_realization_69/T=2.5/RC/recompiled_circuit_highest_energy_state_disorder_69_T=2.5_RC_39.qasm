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
rz(1.5975098) q[0];
sx q[0];
rz(-1.37356) q[0];
sx q[0];
rz(-2.1231667) q[0];
rz(-0.51813689) q[1];
sx q[1];
rz(-2.3845446) q[1];
sx q[1];
rz(0.63176027) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.951183) q[0];
sx q[0];
rz(-2.8993239) q[0];
sx q[0];
rz(-2.1184068) q[0];
rz(-pi) q[1];
rz(1.873513) q[2];
sx q[2];
rz(-1.5945425) q[2];
sx q[2];
rz(2.077092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5093928) q[1];
sx q[1];
rz(-0.94974564) q[1];
sx q[1];
rz(-2.7953447) q[1];
rz(-pi) q[2];
rz(-1.8381836) q[3];
sx q[3];
rz(-2.8041556) q[3];
sx q[3];
rz(-0.22663675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53039256) q[2];
sx q[2];
rz(-0.35553122) q[2];
sx q[2];
rz(1.6543039) q[2];
rz(-0.90855956) q[3];
sx q[3];
rz(-0.23659758) q[3];
sx q[3];
rz(-2.1366185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1210043) q[0];
sx q[0];
rz(-1.4326743) q[0];
sx q[0];
rz(0.16227907) q[0];
rz(0.24457112) q[1];
sx q[1];
rz(-1.9385612) q[1];
sx q[1];
rz(-2.8499106) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2431716) q[0];
sx q[0];
rz(-0.34472877) q[0];
sx q[0];
rz(2.5401462) q[0];
rz(-pi) q[1];
rz(1.2873093) q[2];
sx q[2];
rz(-1.3559196) q[2];
sx q[2];
rz(-0.3706929) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4710992) q[1];
sx q[1];
rz(-2.3355977) q[1];
sx q[1];
rz(-0.32260311) q[1];
rz(-2.347885) q[3];
sx q[3];
rz(-0.55665239) q[3];
sx q[3];
rz(1.8451549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67819277) q[2];
sx q[2];
rz(-2.2726111) q[2];
sx q[2];
rz(-1.0758859) q[2];
rz(-1.894527) q[3];
sx q[3];
rz(-1.5970634) q[3];
sx q[3];
rz(-1.4472848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70188824) q[0];
sx q[0];
rz(-2.0151558) q[0];
sx q[0];
rz(2.9578748) q[0];
rz(1.6366929) q[1];
sx q[1];
rz(-2.3972062) q[1];
sx q[1];
rz(-2.9768129) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713584) q[0];
sx q[0];
rz(-2.5557098) q[0];
sx q[0];
rz(-1.5878994) q[0];
x q[1];
rz(-0.67586502) q[2];
sx q[2];
rz(-1.8177351) q[2];
sx q[2];
rz(0.20154146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3850194) q[1];
sx q[1];
rz(-2.1655271) q[1];
sx q[1];
rz(-1.1043332) q[1];
rz(2.1042473) q[3];
sx q[3];
rz(-1.8672393) q[3];
sx q[3];
rz(-2.8840349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.034721) q[2];
sx q[2];
rz(-1.9811337) q[2];
sx q[2];
rz(2.0351694) q[2];
rz(2.4404081) q[3];
sx q[3];
rz(-1.3283575) q[3];
sx q[3];
rz(-0.4655233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.768854) q[0];
sx q[0];
rz(-2.0059858) q[0];
sx q[0];
rz(-0.94451529) q[0];
rz(1.6230029) q[1];
sx q[1];
rz(-0.98097643) q[1];
sx q[1];
rz(1.6015582) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1683274) q[0];
sx q[0];
rz(-1.0078609) q[0];
sx q[0];
rz(-0.8416147) q[0];
rz(-pi) q[1];
rz(1.4532907) q[2];
sx q[2];
rz(-1.4521993) q[2];
sx q[2];
rz(-0.30726739) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.51241261) q[1];
sx q[1];
rz(-0.80958074) q[1];
sx q[1];
rz(1.045524) q[1];
x q[2];
rz(0.75042689) q[3];
sx q[3];
rz(-1.5380757) q[3];
sx q[3];
rz(1.0270384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1314142) q[2];
sx q[2];
rz(-2.5203036) q[2];
sx q[2];
rz(2.9739001) q[2];
rz(3.0883279) q[3];
sx q[3];
rz(-2.0361418) q[3];
sx q[3];
rz(-1.7648599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5662956) q[0];
sx q[0];
rz(-3.0360041) q[0];
sx q[0];
rz(0.29933023) q[0];
rz(-2.7686367) q[1];
sx q[1];
rz(-1.2495709) q[1];
sx q[1];
rz(-1.6711055) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0051992) q[0];
sx q[0];
rz(-2.1182501) q[0];
sx q[0];
rz(1.0639079) q[0];
rz(-pi) q[1];
rz(-1.6909356) q[2];
sx q[2];
rz(-2.0317626) q[2];
sx q[2];
rz(-1.1663988) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.336323) q[1];
sx q[1];
rz(-2.0741558) q[1];
sx q[1];
rz(0.663228) q[1];
rz(-pi) q[2];
rz(1.3277131) q[3];
sx q[3];
rz(-0.63589261) q[3];
sx q[3];
rz(-2.2896374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2028929) q[2];
sx q[2];
rz(-0.6627658) q[2];
sx q[2];
rz(1.1104442) q[2];
rz(0.86197305) q[3];
sx q[3];
rz(-1.5098666) q[3];
sx q[3];
rz(-1.2601669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2170169) q[0];
sx q[0];
rz(-0.68704263) q[0];
sx q[0];
rz(-2.7365015) q[0];
rz(-0.336126) q[1];
sx q[1];
rz(-1.4210217) q[1];
sx q[1];
rz(-0.95132557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66621214) q[0];
sx q[0];
rz(-2.6586091) q[0];
sx q[0];
rz(-2.5673812) q[0];
rz(-pi) q[1];
rz(-0.73619618) q[2];
sx q[2];
rz(-1.7288107) q[2];
sx q[2];
rz(2.406213) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8321633) q[1];
sx q[1];
rz(-2.1346666) q[1];
sx q[1];
rz(1.840074) q[1];
rz(-1.7634994) q[3];
sx q[3];
rz(-1.8451705) q[3];
sx q[3];
rz(1.5128795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.090791) q[2];
sx q[2];
rz(-1.4046706) q[2];
sx q[2];
rz(0.23078272) q[2];
rz(2.9595621) q[3];
sx q[3];
rz(-2.5735276) q[3];
sx q[3];
rz(2.6128984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7886605) q[0];
sx q[0];
rz(-0.039529888) q[0];
sx q[0];
rz(-0.93210644) q[0];
rz(0.63356361) q[1];
sx q[1];
rz(-1.1195868) q[1];
sx q[1];
rz(1.8036141) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0377696) q[0];
sx q[0];
rz(-1.4892206) q[0];
sx q[0];
rz(1.6768811) q[0];
x q[1];
rz(0.39057486) q[2];
sx q[2];
rz(-2.8468067) q[2];
sx q[2];
rz(0.24402555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89796126) q[1];
sx q[1];
rz(-2.1719031) q[1];
sx q[1];
rz(-2.6578085) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6026434) q[3];
sx q[3];
rz(-1.0726811) q[3];
sx q[3];
rz(-0.3146047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46447095) q[2];
sx q[2];
rz(-1.2268343) q[2];
sx q[2];
rz(-1.9390437) q[2];
rz(2.7770212) q[3];
sx q[3];
rz(-2.4489844) q[3];
sx q[3];
rz(-2.306849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6637591) q[0];
sx q[0];
rz(-0.57357016) q[0];
sx q[0];
rz(0.17663503) q[0];
rz(2.7015576) q[1];
sx q[1];
rz(-1.1853848) q[1];
sx q[1];
rz(0.96493351) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7149219) q[0];
sx q[0];
rz(-1.7650801) q[0];
sx q[0];
rz(1.705709) q[0];
rz(-2.1148483) q[2];
sx q[2];
rz(-2.1026582) q[2];
sx q[2];
rz(0.085467664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31613708) q[1];
sx q[1];
rz(-1.2082005) q[1];
sx q[1];
rz(-0.98939244) q[1];
x q[2];
rz(-0.075841622) q[3];
sx q[3];
rz(-1.8245398) q[3];
sx q[3];
rz(-0.50204078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9289916) q[2];
sx q[2];
rz(-1.5233728) q[2];
sx q[2];
rz(3.1316481) q[2];
rz(-3.0250004) q[3];
sx q[3];
rz(-2.8023585) q[3];
sx q[3];
rz(2.5998083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5803439) q[0];
sx q[0];
rz(-1.2224226) q[0];
sx q[0];
rz(0.62193459) q[0];
rz(-1.4382582) q[1];
sx q[1];
rz(-2.56918) q[1];
sx q[1];
rz(2.4780746) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91231649) q[0];
sx q[0];
rz(-1.7890507) q[0];
sx q[0];
rz(-2.403232) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8542669) q[2];
sx q[2];
rz(-1.4525692) q[2];
sx q[2];
rz(-0.17519874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2199557) q[1];
sx q[1];
rz(-1.8884648) q[1];
sx q[1];
rz(2.2457473) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8264753) q[3];
sx q[3];
rz(-2.8459918) q[3];
sx q[3];
rz(0.35741266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5180987) q[2];
sx q[2];
rz(-1.3891209) q[2];
sx q[2];
rz(0.36250472) q[2];
rz(2.6643961) q[3];
sx q[3];
rz(-2.0878744) q[3];
sx q[3];
rz(1.1227192) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0839194) q[0];
sx q[0];
rz(-2.4102983) q[0];
sx q[0];
rz(0.90173632) q[0];
rz(2.4665191) q[1];
sx q[1];
rz(-2.156064) q[1];
sx q[1];
rz(-2.9973082) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77567277) q[0];
sx q[0];
rz(-2.7808041) q[0];
sx q[0];
rz(-1.3350639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4298019) q[2];
sx q[2];
rz(-0.24082213) q[2];
sx q[2];
rz(-2.4483829) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.044202494) q[1];
sx q[1];
rz(-0.69141885) q[1];
sx q[1];
rz(-0.39908646) q[1];
x q[2];
rz(0.4422632) q[3];
sx q[3];
rz(-1.398842) q[3];
sx q[3];
rz(-2.8124335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5054063) q[2];
sx q[2];
rz(-1.212684) q[2];
sx q[2];
rz(0.20067659) q[2];
rz(1.1096654) q[3];
sx q[3];
rz(-2.6789013) q[3];
sx q[3];
rz(-0.5717352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5889482) q[0];
sx q[0];
rz(-2.1727967) q[0];
sx q[0];
rz(-1.1217242) q[0];
rz(-2.6523392) q[1];
sx q[1];
rz(-1.4514634) q[1];
sx q[1];
rz(-1.0101752) q[1];
rz(-1.1560925) q[2];
sx q[2];
rz(-1.6558052) q[2];
sx q[2];
rz(2.2075352) q[2];
rz(-0.060754178) q[3];
sx q[3];
rz(-2.698632) q[3];
sx q[3];
rz(-2.5198577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
