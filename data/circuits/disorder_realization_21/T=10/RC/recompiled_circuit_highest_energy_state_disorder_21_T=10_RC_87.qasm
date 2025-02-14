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
rz(0.30581185) q[0];
sx q[0];
rz(-2.5078791) q[0];
sx q[0];
rz(0.60854882) q[0];
rz(2.4701056) q[1];
sx q[1];
rz(-0.66317135) q[1];
sx q[1];
rz(1.6821678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072351) q[0];
sx q[0];
rz(-0.57091129) q[0];
sx q[0];
rz(1.1263322) q[0];
rz(-pi) q[1];
rz(1.1205961) q[2];
sx q[2];
rz(-1.3687689) q[2];
sx q[2];
rz(-1.1062619) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8512177) q[1];
sx q[1];
rz(-2.3688815) q[1];
sx q[1];
rz(0.98224838) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3988092) q[3];
sx q[3];
rz(-0.81828749) q[3];
sx q[3];
rz(-2.9060272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7451611) q[2];
sx q[2];
rz(-0.67239422) q[2];
sx q[2];
rz(-0.82610899) q[2];
rz(0.36469001) q[3];
sx q[3];
rz(-1.6148022) q[3];
sx q[3];
rz(2.8472692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.9593338) q[0];
sx q[0];
rz(-1.8987645) q[0];
sx q[0];
rz(-0.44573319) q[0];
rz(2.5511197) q[1];
sx q[1];
rz(-1.0346552) q[1];
sx q[1];
rz(2.7139434) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2954461) q[0];
sx q[0];
rz(-1.5238191) q[0];
sx q[0];
rz(-2.0278005) q[0];
rz(-2.1600464) q[2];
sx q[2];
rz(-2.2950811) q[2];
sx q[2];
rz(-1.8777478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.97274894) q[1];
sx q[1];
rz(-2.4764937) q[1];
sx q[1];
rz(-1.14023) q[1];
x q[2];
rz(-1.5758954) q[3];
sx q[3];
rz(-2.3672464) q[3];
sx q[3];
rz(2.4536595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6815971) q[2];
sx q[2];
rz(-0.57969105) q[2];
sx q[2];
rz(-1.0350636) q[2];
rz(1.5929818) q[3];
sx q[3];
rz(-2.9620453) q[3];
sx q[3];
rz(-2.2523994) q[3];
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
rz(0.11397938) q[0];
sx q[0];
rz(-0.005682156) q[0];
sx q[0];
rz(3.0300544) q[0];
rz(-1.5533718) q[1];
sx q[1];
rz(-0.40451834) q[1];
sx q[1];
rz(-2.9389971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3521393) q[0];
sx q[0];
rz(-1.7048262) q[0];
sx q[0];
rz(1.171597) q[0];
rz(-pi) q[1];
rz(-2.3417046) q[2];
sx q[2];
rz(-2.7383907) q[2];
sx q[2];
rz(-3.0416807) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9249768) q[1];
sx q[1];
rz(-2.1053211) q[1];
sx q[1];
rz(0.85304867) q[1];
x q[2];
rz(-1.9913909) q[3];
sx q[3];
rz(-2.1199825) q[3];
sx q[3];
rz(-2.7830916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69278875) q[2];
sx q[2];
rz(-2.0710129) q[2];
sx q[2];
rz(2.1602574) q[2];
rz(-0.28686178) q[3];
sx q[3];
rz(-0.57932866) q[3];
sx q[3];
rz(2.9441492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.2364748) q[0];
sx q[0];
rz(-3.0906257) q[0];
sx q[0];
rz(-0.68848759) q[0];
rz(0.14432898) q[1];
sx q[1];
rz(-1.1878443) q[1];
sx q[1];
rz(-2.3932638) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9084602) q[0];
sx q[0];
rz(-1.615342) q[0];
sx q[0];
rz(2.4219803) q[0];
rz(-pi) q[1];
rz(-1.8725206) q[2];
sx q[2];
rz(-1.8278484) q[2];
sx q[2];
rz(2.1268877) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1074368) q[1];
sx q[1];
rz(-1.7567051) q[1];
sx q[1];
rz(-0.30349813) q[1];
x q[2];
rz(-0.87012989) q[3];
sx q[3];
rz(-1.8091473) q[3];
sx q[3];
rz(0.42422653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0904842) q[2];
sx q[2];
rz(-1.313831) q[2];
sx q[2];
rz(-1.2468106) q[2];
rz(3.0221853) q[3];
sx q[3];
rz(-1.9236919) q[3];
sx q[3];
rz(-0.30521211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7684105) q[0];
sx q[0];
rz(-0.38653448) q[0];
sx q[0];
rz(-0.61491948) q[0];
rz(-0.83456314) q[1];
sx q[1];
rz(-1.7465697) q[1];
sx q[1];
rz(-0.83647299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1008719) q[0];
sx q[0];
rz(-3.1060069) q[0];
sx q[0];
rz(1.3077223) q[0];
rz(-pi) q[1];
rz(-2.4397298) q[2];
sx q[2];
rz(-2.3387675) q[2];
sx q[2];
rz(-0.35477625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1783318) q[1];
sx q[1];
rz(-1.6815917) q[1];
sx q[1];
rz(-2.1352467) q[1];
x q[2];
rz(2.3079268) q[3];
sx q[3];
rz(-2.1676237) q[3];
sx q[3];
rz(2.4868929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5323083) q[2];
sx q[2];
rz(-2.534635) q[2];
sx q[2];
rz(-0.92030805) q[2];
rz(-2.3406384) q[3];
sx q[3];
rz(-0.52217537) q[3];
sx q[3];
rz(0.33867684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9698708) q[0];
sx q[0];
rz(-1.0665749) q[0];
sx q[0];
rz(-0.26611662) q[0];
rz(-1.7099821) q[1];
sx q[1];
rz(-2.0224729) q[1];
sx q[1];
rz(-2.5439579) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5567096) q[0];
sx q[0];
rz(-2.1765255) q[0];
sx q[0];
rz(-1.8023876) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9605537) q[2];
sx q[2];
rz(-1.5818051) q[2];
sx q[2];
rz(-2.371108) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1823765) q[1];
sx q[1];
rz(-0.87494779) q[1];
sx q[1];
rz(-1.6435502) q[1];
x q[2];
rz(2.9868578) q[3];
sx q[3];
rz(-1.2535963) q[3];
sx q[3];
rz(-1.6586608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0133682) q[2];
sx q[2];
rz(-0.2672264) q[2];
sx q[2];
rz(2.5426148) q[2];
rz(2.525575) q[3];
sx q[3];
rz(-2.7003761) q[3];
sx q[3];
rz(-2.0986572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1004341) q[0];
sx q[0];
rz(-2.2293595) q[0];
sx q[0];
rz(0.31901932) q[0];
rz(3.0306385) q[1];
sx q[1];
rz(-1.3919421) q[1];
sx q[1];
rz(2.9999733) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6017766) q[0];
sx q[0];
rz(-0.68968533) q[0];
sx q[0];
rz(-0.15183555) q[0];
x q[1];
rz(2.1221815) q[2];
sx q[2];
rz(-1.5745609) q[2];
sx q[2];
rz(2.5906445) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95699043) q[1];
sx q[1];
rz(-2.5612368) q[1];
sx q[1];
rz(-0.34837153) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9913031) q[3];
sx q[3];
rz(-1.0553011) q[3];
sx q[3];
rz(-0.084160683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5854599) q[2];
sx q[2];
rz(-1.0254878) q[2];
sx q[2];
rz(0.11036135) q[2];
rz(-0.50802556) q[3];
sx q[3];
rz(-0.042162687) q[3];
sx q[3];
rz(1.149811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13462774) q[0];
sx q[0];
rz(-0.65161324) q[0];
sx q[0];
rz(1.2192669) q[0];
rz(0.02267516) q[1];
sx q[1];
rz(-0.89212787) q[1];
sx q[1];
rz(2.5882744) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0344201) q[0];
sx q[0];
rz(-1.4822472) q[0];
sx q[0];
rz(-2.2293985) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1453317) q[2];
sx q[2];
rz(-1.2424412) q[2];
sx q[2];
rz(2.7789214) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7003498) q[1];
sx q[1];
rz(-1.2257595) q[1];
sx q[1];
rz(-0.54919589) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24711547) q[3];
sx q[3];
rz(-1.860642) q[3];
sx q[3];
rz(-2.4245976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7360709) q[2];
sx q[2];
rz(-0.61554474) q[2];
sx q[2];
rz(1.3496529) q[2];
rz(3.0810629) q[3];
sx q[3];
rz(-0.90977257) q[3];
sx q[3];
rz(0.6282261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26199207) q[0];
sx q[0];
rz(-0.99402004) q[0];
sx q[0];
rz(0.0066198786) q[0];
rz(0.63893843) q[1];
sx q[1];
rz(-0.21955755) q[1];
sx q[1];
rz(2.7117597) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.854786) q[0];
sx q[0];
rz(-1.2116951) q[0];
sx q[0];
rz(2.7451928) q[0];
rz(-pi) q[1];
rz(1.1388135) q[2];
sx q[2];
rz(-1.3393992) q[2];
sx q[2];
rz(-2.9803986) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62334261) q[1];
sx q[1];
rz(-2.3584705) q[1];
sx q[1];
rz(0.030695774) q[1];
x q[2];
rz(0.10363522) q[3];
sx q[3];
rz(-2.5536827) q[3];
sx q[3];
rz(-2.8757446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85475737) q[2];
sx q[2];
rz(-1.1755875) q[2];
sx q[2];
rz(3.1268934) q[2];
rz(1.1356575) q[3];
sx q[3];
rz(-2.0292) q[3];
sx q[3];
rz(1.8358102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056674615) q[0];
sx q[0];
rz(-2.882353) q[0];
sx q[0];
rz(-1.2925451) q[0];
rz(2.9855285) q[1];
sx q[1];
rz(-0.81764692) q[1];
sx q[1];
rz(-2.3110716) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1908295) q[0];
sx q[0];
rz(-2.2939993) q[0];
sx q[0];
rz(-0.13287414) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8409675) q[2];
sx q[2];
rz(-0.81659277) q[2];
sx q[2];
rz(-0.40378324) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76354181) q[1];
sx q[1];
rz(-1.4908067) q[1];
sx q[1];
rz(-1.3536118) q[1];
rz(2.6220406) q[3];
sx q[3];
rz(-0.41659912) q[3];
sx q[3];
rz(-2.8263457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0326651) q[2];
sx q[2];
rz(-0.69643164) q[2];
sx q[2];
rz(-0.44309524) q[2];
rz(2.9923934) q[3];
sx q[3];
rz(-0.34330338) q[3];
sx q[3];
rz(-2.8086015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14015848) q[0];
sx q[0];
rz(-0.59068155) q[0];
sx q[0];
rz(-0.97472192) q[0];
rz(0.99826605) q[1];
sx q[1];
rz(-1.9010192) q[1];
sx q[1];
rz(-0.98767282) q[1];
rz(-0.34299359) q[2];
sx q[2];
rz(-1.7609114) q[2];
sx q[2];
rz(-1.5821725) q[2];
rz(2.0988437) q[3];
sx q[3];
rz(-0.73560148) q[3];
sx q[3];
rz(-0.92675496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
