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
rz(-0.15859088) q[0];
sx q[0];
rz(2.6245485) q[0];
sx q[0];
rz(11.256097) q[0];
rz(-0.34215555) q[1];
sx q[1];
rz(-1.181239) q[1];
sx q[1];
rz(-2.1769843) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83822891) q[0];
sx q[0];
rz(-1.705929) q[0];
sx q[0];
rz(-2.5101296) q[0];
x q[1];
rz(2.7959441) q[2];
sx q[2];
rz(-1.3866349) q[2];
sx q[2];
rz(-0.08376567) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55935153) q[1];
sx q[1];
rz(-0.808945) q[1];
sx q[1];
rz(-1.9870583) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11537376) q[3];
sx q[3];
rz(-2.2386207) q[3];
sx q[3];
rz(1.3504207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0059119314) q[2];
sx q[2];
rz(-2.9475806) q[2];
sx q[2];
rz(1.4646336) q[2];
rz(-1.3095193) q[3];
sx q[3];
rz(-0.94683164) q[3];
sx q[3];
rz(2.8215698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2784073) q[0];
sx q[0];
rz(-1.139737) q[0];
sx q[0];
rz(-1.3673258) q[0];
rz(-1.4890081) q[1];
sx q[1];
rz(-2.3431578) q[1];
sx q[1];
rz(0.29261011) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.145702) q[0];
sx q[0];
rz(-1.2371024) q[0];
sx q[0];
rz(-1.6745553) q[0];
rz(-0.16186773) q[2];
sx q[2];
rz(-1.4449931) q[2];
sx q[2];
rz(2.6243072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.039174883) q[1];
sx q[1];
rz(-0.31512773) q[1];
sx q[1];
rz(2.4246895) q[1];
rz(-pi) q[2];
rz(-0.23789101) q[3];
sx q[3];
rz(-2.8679016) q[3];
sx q[3];
rz(-2.9546933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1613529) q[2];
sx q[2];
rz(-2.5148401) q[2];
sx q[2];
rz(-0.18388595) q[2];
rz(-2.6889177) q[3];
sx q[3];
rz(-1.5225182) q[3];
sx q[3];
rz(-3.0116853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9420796) q[0];
sx q[0];
rz(-0.71490723) q[0];
sx q[0];
rz(-2.3364501) q[0];
rz(-1.0792271) q[1];
sx q[1];
rz(-2.3715623) q[1];
sx q[1];
rz(1.315518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41322177) q[0];
sx q[0];
rz(-2.411619) q[0];
sx q[0];
rz(3.0527924) q[0];
rz(-1.075885) q[2];
sx q[2];
rz(-0.80782164) q[2];
sx q[2];
rz(-0.13208315) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6650369) q[1];
sx q[1];
rz(-1.8522693) q[1];
sx q[1];
rz(3.046077) q[1];
x q[2];
rz(1.6571548) q[3];
sx q[3];
rz(-1.8724073) q[3];
sx q[3];
rz(-0.82179606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9341854) q[2];
sx q[2];
rz(-1.4855874) q[2];
sx q[2];
rz(2.5918813) q[2];
rz(-3.1304729) q[3];
sx q[3];
rz(-0.79922262) q[3];
sx q[3];
rz(-1.3219249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84871197) q[0];
sx q[0];
rz(-0.95893812) q[0];
sx q[0];
rz(0.30353656) q[0];
rz(-0.15829463) q[1];
sx q[1];
rz(-2.0469432) q[1];
sx q[1];
rz(2.1319481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15561129) q[0];
sx q[0];
rz(-1.2878875) q[0];
sx q[0];
rz(-2.020218) q[0];
rz(-pi) q[1];
rz(2.7042806) q[2];
sx q[2];
rz(-0.6491937) q[2];
sx q[2];
rz(-1.5199666) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4588488) q[1];
sx q[1];
rz(-1.7284596) q[1];
sx q[1];
rz(2.661652) q[1];
rz(-pi) q[2];
rz(1.5916274) q[3];
sx q[3];
rz(-2.0685141) q[3];
sx q[3];
rz(1.8128943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.9369072) q[2];
sx q[2];
rz(-2.3202809) q[2];
sx q[2];
rz(2.8724907) q[2];
rz(-1.6875632) q[3];
sx q[3];
rz(-1.5917835) q[3];
sx q[3];
rz(-1.5788186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27075416) q[0];
sx q[0];
rz(-2.0814867) q[0];
sx q[0];
rz(-0.36218542) q[0];
rz(-0.65131342) q[1];
sx q[1];
rz(-1.6505227) q[1];
sx q[1];
rz(-1.5247033) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7633879) q[0];
sx q[0];
rz(-1.4414865) q[0];
sx q[0];
rz(0.61152258) q[0];
rz(-pi) q[1];
rz(3.0595018) q[2];
sx q[2];
rz(-2.1920648) q[2];
sx q[2];
rz(2.4022849) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0514879) q[1];
sx q[1];
rz(-1.0516775) q[1];
sx q[1];
rz(0.41579065) q[1];
x q[2];
rz(1.9664496) q[3];
sx q[3];
rz(-0.86651245) q[3];
sx q[3];
rz(-0.3053529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0321956) q[2];
sx q[2];
rz(-0.95462644) q[2];
sx q[2];
rz(1.7677914) q[2];
rz(2.7023756) q[3];
sx q[3];
rz(-0.64911157) q[3];
sx q[3];
rz(-0.60905987) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1010901) q[0];
sx q[0];
rz(-2.4093565) q[0];
sx q[0];
rz(2.8771583) q[0];
rz(-0.6791555) q[1];
sx q[1];
rz(-0.87398386) q[1];
sx q[1];
rz(-0.98495475) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.934498) q[0];
sx q[0];
rz(-2.6885928) q[0];
sx q[0];
rz(2.3191602) q[0];
x q[1];
rz(2.856065) q[2];
sx q[2];
rz(-1.1368903) q[2];
sx q[2];
rz(-1.3615695) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.15036525) q[1];
sx q[1];
rz(-2.5179407) q[1];
sx q[1];
rz(0.6536478) q[1];
rz(-pi) q[2];
rz(1.1285176) q[3];
sx q[3];
rz(-1.3849713) q[3];
sx q[3];
rz(3.1215661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1017477) q[2];
sx q[2];
rz(-0.16182772) q[2];
sx q[2];
rz(2.7092773) q[2];
rz(1.013422) q[3];
sx q[3];
rz(-1.6067959) q[3];
sx q[3];
rz(2.5920946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48563114) q[0];
sx q[0];
rz(-0.39325842) q[0];
sx q[0];
rz(-1.1095169) q[0];
rz(-2.3484777) q[1];
sx q[1];
rz(-1.3536645) q[1];
sx q[1];
rz(-1.5951593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6731747) q[0];
sx q[0];
rz(-2.1270848) q[0];
sx q[0];
rz(-0.89222096) q[0];
x q[1];
rz(1.7810473) q[2];
sx q[2];
rz(-2.2814007) q[2];
sx q[2];
rz(-0.17786121) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7132414) q[1];
sx q[1];
rz(-2.4292786) q[1];
sx q[1];
rz(-1.8800432) q[1];
rz(-1.6768544) q[3];
sx q[3];
rz(-2.2844188) q[3];
sx q[3];
rz(-2.3372363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82211295) q[2];
sx q[2];
rz(-0.18222465) q[2];
sx q[2];
rz(-1.7721843) q[2];
rz(-2.2111514) q[3];
sx q[3];
rz(-1.7934099) q[3];
sx q[3];
rz(1.5038917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.9419788) q[0];
sx q[0];
rz(-1.9604585) q[0];
sx q[0];
rz(0.30366316) q[0];
rz(-2.1619201) q[1];
sx q[1];
rz(-0.62398463) q[1];
sx q[1];
rz(-0.40506515) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41202085) q[0];
sx q[0];
rz(-1.3776017) q[0];
sx q[0];
rz(-0.3839107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.151565) q[2];
sx q[2];
rz(-1.4128601) q[2];
sx q[2];
rz(-1.4106223) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77877264) q[1];
sx q[1];
rz(-1.3937794) q[1];
sx q[1];
rz(2.5702012) q[1];
rz(-pi) q[2];
rz(1.1438683) q[3];
sx q[3];
rz(-1.3865304) q[3];
sx q[3];
rz(-1.1729405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27641174) q[2];
sx q[2];
rz(-0.96651912) q[2];
sx q[2];
rz(2.442404) q[2];
rz(0.80896038) q[3];
sx q[3];
rz(-1.304108) q[3];
sx q[3];
rz(2.8185524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9998099) q[0];
sx q[0];
rz(-1.1572105) q[0];
sx q[0];
rz(0.13741563) q[0];
rz(-3.0925062) q[1];
sx q[1];
rz(-2.0259435) q[1];
sx q[1];
rz(-0.5079937) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8592472) q[0];
sx q[0];
rz(-1.0758335) q[0];
sx q[0];
rz(-1.540394) q[0];
rz(-1.5586583) q[2];
sx q[2];
rz(-1.3985139) q[2];
sx q[2];
rz(0.42086312) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0145871) q[1];
sx q[1];
rz(-2.0620605) q[1];
sx q[1];
rz(0.36775695) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6925595) q[3];
sx q[3];
rz(-1.4819488) q[3];
sx q[3];
rz(0.73032398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5892443) q[2];
sx q[2];
rz(-1.9231223) q[2];
sx q[2];
rz(0.98769665) q[2];
rz(-0.47932953) q[3];
sx q[3];
rz(-1.5440116) q[3];
sx q[3];
rz(0.89232579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.669303) q[0];
sx q[0];
rz(-0.23858128) q[0];
sx q[0];
rz(-3.0008089) q[0];
rz(-2.7760778) q[1];
sx q[1];
rz(-0.82217685) q[1];
sx q[1];
rz(2.4929094) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4869944) q[0];
sx q[0];
rz(-2.2077422) q[0];
sx q[0];
rz(1.1331228) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12154433) q[2];
sx q[2];
rz(-0.76833188) q[2];
sx q[2];
rz(-0.8476846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59274793) q[1];
sx q[1];
rz(-0.59863585) q[1];
sx q[1];
rz(3.0416802) q[1];
rz(-2.8148531) q[3];
sx q[3];
rz(-1.3413197) q[3];
sx q[3];
rz(2.6658201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4623798) q[2];
sx q[2];
rz(-2.1361735) q[2];
sx q[2];
rz(0.072754808) q[2];
rz(0.66271979) q[3];
sx q[3];
rz(-1.8279671) q[3];
sx q[3];
rz(-0.16547468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0543906) q[0];
sx q[0];
rz(-1.6086171) q[0];
sx q[0];
rz(1.4679012) q[0];
rz(-1.5351334) q[1];
sx q[1];
rz(-2.2538593) q[1];
sx q[1];
rz(-1.5182553) q[1];
rz(0.4735422) q[2];
sx q[2];
rz(-2.8859856) q[2];
sx q[2];
rz(-2.4398266) q[2];
rz(2.0562511) q[3];
sx q[3];
rz(-2.3997636) q[3];
sx q[3];
rz(2.6482481) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
