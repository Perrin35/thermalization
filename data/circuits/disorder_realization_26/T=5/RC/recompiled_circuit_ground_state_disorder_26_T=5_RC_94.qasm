OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3684664) q[0];
sx q[0];
rz(-2.3946895) q[0];
sx q[0];
rz(-0.84063831) q[0];
rz(0.12159881) q[1];
sx q[1];
rz(-1.2727979) q[1];
sx q[1];
rz(0.29903856) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7515669) q[0];
sx q[0];
rz(-1.9932858) q[0];
sx q[0];
rz(2.3038191) q[0];
rz(-1.0844803) q[2];
sx q[2];
rz(-1.63967) q[2];
sx q[2];
rz(-1.8191847) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2192397) q[1];
sx q[1];
rz(-1.427622) q[1];
sx q[1];
rz(2.8994865) q[1];
rz(1.8567919) q[3];
sx q[3];
rz(-2.6702849) q[3];
sx q[3];
rz(-2.907674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76089871) q[2];
sx q[2];
rz(-2.0913405) q[2];
sx q[2];
rz(-2.8139581) q[2];
rz(-1.7662175) q[3];
sx q[3];
rz(-1.4204493) q[3];
sx q[3];
rz(-0.49427858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7647917) q[0];
sx q[0];
rz(-1.5400274) q[0];
sx q[0];
rz(0.85533992) q[0];
rz(0.0080464706) q[1];
sx q[1];
rz(-1.2866373) q[1];
sx q[1];
rz(-2.6904552) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9066042) q[0];
sx q[0];
rz(-1.7275817) q[0];
sx q[0];
rz(-2.2995728) q[0];
rz(-pi) q[1];
rz(-0.77215423) q[2];
sx q[2];
rz(-1.4322965) q[2];
sx q[2];
rz(2.0318299) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5246213) q[1];
sx q[1];
rz(-2.0071623) q[1];
sx q[1];
rz(1.7625767) q[1];
rz(-pi) q[2];
rz(0.69591095) q[3];
sx q[3];
rz(-0.59248052) q[3];
sx q[3];
rz(-2.185948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1796639) q[2];
sx q[2];
rz(-2.0018061) q[2];
sx q[2];
rz(-3.0120604) q[2];
rz(-2.9642963) q[3];
sx q[3];
rz(-0.57759053) q[3];
sx q[3];
rz(-3.0887443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.391908) q[0];
sx q[0];
rz(-2.7295697) q[0];
sx q[0];
rz(2.5823197) q[0];
rz(0.094712146) q[1];
sx q[1];
rz(-1.5385224) q[1];
sx q[1];
rz(-2.6643378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1806948) q[0];
sx q[0];
rz(-1.8261693) q[0];
sx q[0];
rz(-1.1605422) q[0];
rz(-1.4213218) q[2];
sx q[2];
rz(-1.9286619) q[2];
sx q[2];
rz(-2.0400816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4871769) q[1];
sx q[1];
rz(-0.80779874) q[1];
sx q[1];
rz(0.31262763) q[1];
rz(0.52366728) q[3];
sx q[3];
rz(-0.91991495) q[3];
sx q[3];
rz(1.6802579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7462848) q[2];
sx q[2];
rz(-1.8751273) q[2];
sx q[2];
rz(-0.9300119) q[2];
rz(1.8837455) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(0.045684489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3607218) q[0];
sx q[0];
rz(-1.5683132) q[0];
sx q[0];
rz(-1.038653) q[0];
rz(-0.61141283) q[1];
sx q[1];
rz(-0.88714209) q[1];
sx q[1];
rz(0.92322737) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99902376) q[0];
sx q[0];
rz(-1.9143189) q[0];
sx q[0];
rz(-0.96238636) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.011767894) q[2];
sx q[2];
rz(-1.1530877) q[2];
sx q[2];
rz(-2.2209446) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23455305) q[1];
sx q[1];
rz(-0.72826339) q[1];
sx q[1];
rz(-1.6978463) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3096894) q[3];
sx q[3];
rz(-0.93702836) q[3];
sx q[3];
rz(1.5719617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5103147) q[2];
sx q[2];
rz(-2.6738561) q[2];
sx q[2];
rz(2.5941217) q[2];
rz(-0.064420961) q[3];
sx q[3];
rz(-1.3037325) q[3];
sx q[3];
rz(-2.2678383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9861458) q[0];
sx q[0];
rz(-0.90279818) q[0];
sx q[0];
rz(1.4601532) q[0];
rz(-2.1414781) q[1];
sx q[1];
rz(-0.8668879) q[1];
sx q[1];
rz(3.0217081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45792199) q[0];
sx q[0];
rz(-0.73212762) q[0];
sx q[0];
rz(-1.9089655) q[0];
x q[1];
rz(0.5849383) q[2];
sx q[2];
rz(-2.1368933) q[2];
sx q[2];
rz(1.0370129) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2561349) q[1];
sx q[1];
rz(-1.7856344) q[1];
sx q[1];
rz(1.6515405) q[1];
rz(-pi) q[2];
rz(-0.67976953) q[3];
sx q[3];
rz(-2.418737) q[3];
sx q[3];
rz(0.86833176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.039006058) q[2];
sx q[2];
rz(-2.234499) q[2];
sx q[2];
rz(0.26068035) q[2];
rz(-2.2634704) q[3];
sx q[3];
rz(-1.2295281) q[3];
sx q[3];
rz(-2.3584283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6108625) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(2.1395785) q[0];
rz(-0.37725457) q[1];
sx q[1];
rz(-0.95163029) q[1];
sx q[1];
rz(2.196905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2827891) q[0];
sx q[0];
rz(-1.3529881) q[0];
sx q[0];
rz(-2.5522638) q[0];
x q[1];
rz(-2.2425507) q[2];
sx q[2];
rz(-2.0625249) q[2];
sx q[2];
rz(1.1373625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73026555) q[1];
sx q[1];
rz(-1.7217727) q[1];
sx q[1];
rz(-0.32457268) q[1];
rz(1.3768436) q[3];
sx q[3];
rz(-2.5911281) q[3];
sx q[3];
rz(2.6725519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79823309) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(-0.043924335) q[2];
rz(-1.6262866) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(1.4656434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8632904) q[0];
sx q[0];
rz(-0.55178061) q[0];
sx q[0];
rz(2.5569051) q[0];
rz(2.0901285) q[1];
sx q[1];
rz(-0.81948558) q[1];
sx q[1];
rz(-0.47031602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3669339) q[0];
sx q[0];
rz(-0.99307382) q[0];
sx q[0];
rz(-1.7455186) q[0];
rz(-2.7696848) q[2];
sx q[2];
rz(-1.1924517) q[2];
sx q[2];
rz(-1.3828948) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9661588) q[1];
sx q[1];
rz(-2.325238) q[1];
sx q[1];
rz(2.6360126) q[1];
rz(-pi) q[2];
rz(-3.028585) q[3];
sx q[3];
rz(-1.5829493) q[3];
sx q[3];
rz(3.0719824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26479244) q[2];
sx q[2];
rz(-2.6046643) q[2];
sx q[2];
rz(2.8301767) q[2];
rz(0.16820678) q[3];
sx q[3];
rz(-1.5663389) q[3];
sx q[3];
rz(0.28347191) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6300221) q[0];
sx q[0];
rz(-0.011307414) q[0];
sx q[0];
rz(-2.2054963) q[0];
rz(-2.6452737) q[1];
sx q[1];
rz(-2.4744108) q[1];
sx q[1];
rz(2.671303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59018001) q[0];
sx q[0];
rz(-1.0820933) q[0];
sx q[0];
rz(-1.2312698) q[0];
x q[1];
rz(3.1390433) q[2];
sx q[2];
rz(-1.7685206) q[2];
sx q[2];
rz(1.6127123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0691904) q[1];
sx q[1];
rz(-1.6399929) q[1];
sx q[1];
rz(1.2284159) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7467198) q[3];
sx q[3];
rz(-1.298549) q[3];
sx q[3];
rz(1.7737089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5101667) q[2];
sx q[2];
rz(-1.7743856) q[2];
sx q[2];
rz(-1.6660956) q[2];
rz(-0.79536074) q[3];
sx q[3];
rz(-0.16203351) q[3];
sx q[3];
rz(1.2696666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0610166) q[0];
sx q[0];
rz(-2.5503655) q[0];
sx q[0];
rz(-3.0812145) q[0];
rz(-2.9810442) q[1];
sx q[1];
rz(-1.5269273) q[1];
sx q[1];
rz(-2.1626332) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706021) q[0];
sx q[0];
rz(-2.5033062) q[0];
sx q[0];
rz(2.9459475) q[0];
x q[1];
rz(0.27411119) q[2];
sx q[2];
rz(-2.2345671) q[2];
sx q[2];
rz(1.2520777) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95293249) q[1];
sx q[1];
rz(-0.97129909) q[1];
sx q[1];
rz(-0.99296928) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5185202) q[3];
sx q[3];
rz(-0.57988047) q[3];
sx q[3];
rz(-1.9573905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0742566) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(-0.072889797) q[2];
rz(0.59761754) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(-0.0028751956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.446796) q[0];
sx q[0];
rz(-2.1294761) q[0];
sx q[0];
rz(-2.6126675) q[0];
rz(0.18813285) q[1];
sx q[1];
rz(-2.4362322) q[1];
sx q[1];
rz(2.5573152) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0952632) q[0];
sx q[0];
rz(-1.8285311) q[0];
sx q[0];
rz(-1.3168174) q[0];
rz(0.71888098) q[2];
sx q[2];
rz(-2.5846057) q[2];
sx q[2];
rz(0.62918909) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1297785) q[1];
sx q[1];
rz(-0.41167694) q[1];
sx q[1];
rz(-3.0342388) q[1];
rz(-pi) q[2];
rz(2.1891865) q[3];
sx q[3];
rz(-1.795553) q[3];
sx q[3];
rz(0.99638961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75954413) q[2];
sx q[2];
rz(-1.0039165) q[2];
sx q[2];
rz(-3.0697401) q[2];
rz(2.1182649) q[3];
sx q[3];
rz(-1.5691248) q[3];
sx q[3];
rz(0.48880997) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1857984) q[0];
sx q[0];
rz(-0.6577984) q[0];
sx q[0];
rz(-0.26554769) q[0];
rz(0.67509782) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(2.6266392) q[2];
sx q[2];
rz(-1.6390159) q[2];
sx q[2];
rz(-1.9255571) q[2];
rz(-1.731338) q[3];
sx q[3];
rz(-1.7358801) q[3];
sx q[3];
rz(0.87270234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
