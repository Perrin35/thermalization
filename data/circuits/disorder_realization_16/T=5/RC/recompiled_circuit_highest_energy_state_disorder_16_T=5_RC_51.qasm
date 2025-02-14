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
rz(1.3311812) q[0];
sx q[0];
rz(-1.5067195) q[0];
sx q[0];
rz(-2.9659086) q[0];
rz(-2.076258) q[1];
sx q[1];
rz(-1.330436) q[1];
sx q[1];
rz(0.65520823) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5979008) q[0];
sx q[0];
rz(-0.59698623) q[0];
sx q[0];
rz(2.1587121) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51061012) q[2];
sx q[2];
rz(-3.0282927) q[2];
sx q[2];
rz(2.4093546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44819278) q[1];
sx q[1];
rz(-1.498135) q[1];
sx q[1];
rz(3.0334515) q[1];
x q[2];
rz(-1.9780065) q[3];
sx q[3];
rz(-2.4904479) q[3];
sx q[3];
rz(-2.4442087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6931927) q[2];
sx q[2];
rz(-0.11036631) q[2];
sx q[2];
rz(-2.8006862) q[2];
rz(2.36813) q[3];
sx q[3];
rz(-1.0186467) q[3];
sx q[3];
rz(-3.036518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.676749) q[0];
sx q[0];
rz(-2.8571547) q[0];
sx q[0];
rz(-0.25962096) q[0];
rz(-0.78850857) q[1];
sx q[1];
rz(-2.2854243) q[1];
sx q[1];
rz(1.8964918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7126078) q[0];
sx q[0];
rz(-0.95604205) q[0];
sx q[0];
rz(0.12354596) q[0];
rz(-1.670426) q[2];
sx q[2];
rz(-1.2567695) q[2];
sx q[2];
rz(-3.0977566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4739405) q[1];
sx q[1];
rz(-2.2721935) q[1];
sx q[1];
rz(-0.46469537) q[1];
x q[2];
rz(-2.63381) q[3];
sx q[3];
rz(-1.4243817) q[3];
sx q[3];
rz(-3.0203569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9022687) q[2];
sx q[2];
rz(-2.6332492) q[2];
sx q[2];
rz(2.2321205) q[2];
rz(-2.4746312) q[3];
sx q[3];
rz(-1.4374377) q[3];
sx q[3];
rz(-0.10233574) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42720902) q[0];
sx q[0];
rz(-0.41218555) q[0];
sx q[0];
rz(1.7538196) q[0];
rz(-2.1448403) q[1];
sx q[1];
rz(-2.4501188) q[1];
sx q[1];
rz(2.6549285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5202058) q[0];
sx q[0];
rz(-0.4356111) q[0];
sx q[0];
rz(1.9389047) q[0];
rz(-0.83618323) q[2];
sx q[2];
rz(-2.0942405) q[2];
sx q[2];
rz(1.5627031) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0906252) q[1];
sx q[1];
rz(-2.087346) q[1];
sx q[1];
rz(-1.4314639) q[1];
x q[2];
rz(-2.1975945) q[3];
sx q[3];
rz(-1.6907915) q[3];
sx q[3];
rz(-1.3809518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3643058) q[2];
sx q[2];
rz(-1.2591668) q[2];
sx q[2];
rz(-0.49059179) q[2];
rz(-2.3250438) q[3];
sx q[3];
rz(-2.3994763) q[3];
sx q[3];
rz(-2.1013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8051324) q[0];
sx q[0];
rz(-1.8959624) q[0];
sx q[0];
rz(0.98139393) q[0];
rz(-2.7384752) q[1];
sx q[1];
rz(-2.6336481) q[1];
sx q[1];
rz(2.6037762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11975461) q[0];
sx q[0];
rz(-1.7203385) q[0];
sx q[0];
rz(-1.3215905) q[0];
x q[1];
rz(0.85015202) q[2];
sx q[2];
rz(-1.6714897) q[2];
sx q[2];
rz(3.0792011) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0672978) q[1];
sx q[1];
rz(-2.5038669) q[1];
sx q[1];
rz(-2.2504351) q[1];
rz(-pi) q[2];
rz(0.31258055) q[3];
sx q[3];
rz(-2.2214032) q[3];
sx q[3];
rz(1.3530817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2526523) q[2];
sx q[2];
rz(-0.50709587) q[2];
sx q[2];
rz(2.0895152) q[2];
rz(2.5893411) q[3];
sx q[3];
rz(-1.735732) q[3];
sx q[3];
rz(1.2205869) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1110765) q[0];
sx q[0];
rz(-1.3330326) q[0];
sx q[0];
rz(0.037574969) q[0];
rz(-2.3887718) q[1];
sx q[1];
rz(-1.5767187) q[1];
sx q[1];
rz(-3.0573696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1128453) q[0];
sx q[0];
rz(-1.0431093) q[0];
sx q[0];
rz(-0.28808388) q[0];
rz(-pi) q[1];
rz(1.8657612) q[2];
sx q[2];
rz(-1.9199445) q[2];
sx q[2];
rz(-2.5087506) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4956717) q[1];
sx q[1];
rz(-2.715647) q[1];
sx q[1];
rz(1.4703712) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4100759) q[3];
sx q[3];
rz(-1.525021) q[3];
sx q[3];
rz(-1.9327527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61421824) q[2];
sx q[2];
rz(-2.1671593) q[2];
sx q[2];
rz(0.43240377) q[2];
rz(-1.6480986) q[3];
sx q[3];
rz(-1.7688388) q[3];
sx q[3];
rz(2.4026332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.768141) q[0];
sx q[0];
rz(-1.5274763) q[0];
sx q[0];
rz(-2.3798808) q[0];
rz(-2.1181882) q[1];
sx q[1];
rz(-2.6506212) q[1];
sx q[1];
rz(0.34058079) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5368333) q[0];
sx q[0];
rz(-0.7927466) q[0];
sx q[0];
rz(-1.6954846) q[0];
x q[1];
rz(1.7271785) q[2];
sx q[2];
rz(-2.1264653) q[2];
sx q[2];
rz(1.8506236) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.611515) q[1];
sx q[1];
rz(-2.7002091) q[1];
sx q[1];
rz(-2.2468021) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4177001) q[3];
sx q[3];
rz(-0.82624895) q[3];
sx q[3];
rz(-1.5797256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1731825) q[2];
sx q[2];
rz(-1.1334271) q[2];
sx q[2];
rz(-0.75562149) q[2];
rz(0.38718811) q[3];
sx q[3];
rz(-1.7500992) q[3];
sx q[3];
rz(-0.79425991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25283915) q[0];
sx q[0];
rz(-0.08789739) q[0];
sx q[0];
rz(-1.2095691) q[0];
rz(-1.5225438) q[1];
sx q[1];
rz(-2.3698898) q[1];
sx q[1];
rz(-1.3873772) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6722058) q[0];
sx q[0];
rz(-1.9090561) q[0];
sx q[0];
rz(-0.32899186) q[0];
rz(-pi) q[1];
rz(2.0389236) q[2];
sx q[2];
rz(-2.515677) q[2];
sx q[2];
rz(0.70410888) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55303326) q[1];
sx q[1];
rz(-0.84399763) q[1];
sx q[1];
rz(-0.56808205) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58297662) q[3];
sx q[3];
rz(-0.54123961) q[3];
sx q[3];
rz(1.6868633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9428955) q[2];
sx q[2];
rz(-0.2936475) q[2];
sx q[2];
rz(-0.4100619) q[2];
rz(-0.53747082) q[3];
sx q[3];
rz(-2.0987857) q[3];
sx q[3];
rz(-1.0679831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-3.033567) q[0];
sx q[0];
rz(-1.2355868) q[0];
sx q[0];
rz(2.1349452) q[0];
rz(2.0058696) q[1];
sx q[1];
rz(-2.6508811) q[1];
sx q[1];
rz(0.27166414) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1170809) q[0];
sx q[0];
rz(-2.4371464) q[0];
sx q[0];
rz(-0.68450494) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92046787) q[2];
sx q[2];
rz(-2.1365385) q[2];
sx q[2];
rz(-0.56220245) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.063273009) q[1];
sx q[1];
rz(-2.1542179) q[1];
sx q[1];
rz(-2.2607445) q[1];
x q[2];
rz(1.7944505) q[3];
sx q[3];
rz(-2.3021379) q[3];
sx q[3];
rz(-0.37089965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2037619) q[2];
sx q[2];
rz(-2.2088642) q[2];
sx q[2];
rz(1.1057828) q[2];
rz(1.743099) q[3];
sx q[3];
rz(-2.1574557) q[3];
sx q[3];
rz(0.85247803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6451013) q[0];
sx q[0];
rz(-0.95106769) q[0];
sx q[0];
rz(-0.35354653) q[0];
rz(1.4191267) q[1];
sx q[1];
rz(-1.5357693) q[1];
sx q[1];
rz(0.062573418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4434011) q[0];
sx q[0];
rz(-2.0581322) q[0];
sx q[0];
rz(2.5145825) q[0];
x q[1];
rz(1.0163505) q[2];
sx q[2];
rz(-2.2502463) q[2];
sx q[2];
rz(2.0959977) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8282916) q[1];
sx q[1];
rz(-2.5646659) q[1];
sx q[1];
rz(0.58842701) q[1];
rz(-pi) q[2];
rz(1.6639641) q[3];
sx q[3];
rz(-0.9103295) q[3];
sx q[3];
rz(-2.4371393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.320437) q[2];
sx q[2];
rz(-1.374036) q[2];
sx q[2];
rz(2.721411) q[2];
rz(-0.34475103) q[3];
sx q[3];
rz(-0.77778608) q[3];
sx q[3];
rz(-2.2043118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6688113) q[0];
sx q[0];
rz(-0.16952276) q[0];
sx q[0];
rz(1.8602759) q[0];
rz(-1.5877113) q[1];
sx q[1];
rz(-0.80931598) q[1];
sx q[1];
rz(-0.70714998) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4879935) q[0];
sx q[0];
rz(-0.63051134) q[0];
sx q[0];
rz(1.5911908) q[0];
rz(0.94790876) q[2];
sx q[2];
rz(-0.68205183) q[2];
sx q[2];
rz(-0.94833224) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.36701074) q[1];
sx q[1];
rz(-0.38644192) q[1];
sx q[1];
rz(-2.5143753) q[1];
rz(-0.78451802) q[3];
sx q[3];
rz(-2.3459179) q[3];
sx q[3];
rz(2.8592542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82004929) q[2];
sx q[2];
rz(-2.0115435) q[2];
sx q[2];
rz(2.6450487) q[2];
rz(-1.8494362) q[3];
sx q[3];
rz(-2.8437331) q[3];
sx q[3];
rz(-2.7636102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15688607) q[0];
sx q[0];
rz(-2.8732185) q[0];
sx q[0];
rz(2.1068841) q[0];
rz(2.5490419) q[1];
sx q[1];
rz(-0.68065803) q[1];
sx q[1];
rz(-1.5417644) q[1];
rz(1.8027923) q[2];
sx q[2];
rz(-1.5314039) q[2];
sx q[2];
rz(1.0614131) q[2];
rz(-1.3453398) q[3];
sx q[3];
rz(-0.97889789) q[3];
sx q[3];
rz(-1.2759664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
