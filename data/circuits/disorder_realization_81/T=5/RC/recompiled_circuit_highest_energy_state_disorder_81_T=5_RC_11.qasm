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
rz(0.45646271) q[0];
sx q[0];
rz(-2.2678092) q[0];
sx q[0];
rz(1.1377347) q[0];
rz(0.046009215) q[1];
sx q[1];
rz(-0.53541056) q[1];
sx q[1];
rz(1.6028264) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7932324) q[0];
sx q[0];
rz(-1.5377511) q[0];
sx q[0];
rz(-1.5263625) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2884772) q[2];
sx q[2];
rz(-1.4783646) q[2];
sx q[2];
rz(1.5866367) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0823209) q[1];
sx q[1];
rz(-0.50283018) q[1];
sx q[1];
rz(0.4362274) q[1];
x q[2];
rz(-0.055786415) q[3];
sx q[3];
rz(-1.7100289) q[3];
sx q[3];
rz(3.1007953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5662745) q[2];
sx q[2];
rz(-0.34276572) q[2];
sx q[2];
rz(-2.9052367) q[2];
rz(0.44860873) q[3];
sx q[3];
rz(-1.6581422) q[3];
sx q[3];
rz(-1.0033222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3818632) q[0];
sx q[0];
rz(-0.48668447) q[0];
sx q[0];
rz(0.78616649) q[0];
rz(-0.78474125) q[1];
sx q[1];
rz(-1.6678383) q[1];
sx q[1];
rz(-3.0395708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1624296) q[0];
sx q[0];
rz(-1.7722881) q[0];
sx q[0];
rz(-0.54690806) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29709105) q[2];
sx q[2];
rz(-0.80657437) q[2];
sx q[2];
rz(-1.7809803) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35013546) q[1];
sx q[1];
rz(-1.268549) q[1];
sx q[1];
rz(-1.8552029) q[1];
rz(-pi) q[2];
rz(0.24436538) q[3];
sx q[3];
rz(-0.8042585) q[3];
sx q[3];
rz(-1.3028631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.027779) q[2];
sx q[2];
rz(-1.4872097) q[2];
sx q[2];
rz(2.8786744) q[2];
rz(-1.3024088) q[3];
sx q[3];
rz(-1.115256) q[3];
sx q[3];
rz(2.3579679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0021492783) q[0];
sx q[0];
rz(-0.0037010598) q[0];
sx q[0];
rz(-1.8078467) q[0];
rz(1.7568781) q[1];
sx q[1];
rz(-2.2127071) q[1];
sx q[1];
rz(2.3768916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36604213) q[0];
sx q[0];
rz(-1.8100272) q[0];
sx q[0];
rz(1.3711098) q[0];
rz(-pi) q[1];
rz(-0.3777067) q[2];
sx q[2];
rz(-1.3481082) q[2];
sx q[2];
rz(-1.4770222) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0955268) q[1];
sx q[1];
rz(-0.30991677) q[1];
sx q[1];
rz(-2.168783) q[1];
x q[2];
rz(1.9165048) q[3];
sx q[3];
rz(-2.3937493) q[3];
sx q[3];
rz(1.9077099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9646405) q[2];
sx q[2];
rz(-1.2531589) q[2];
sx q[2];
rz(0.18356744) q[2];
rz(2.99672) q[3];
sx q[3];
rz(-1.262007) q[3];
sx q[3];
rz(-0.22411331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(2.1862828) q[0];
sx q[0];
rz(-2.0109542) q[0];
sx q[0];
rz(0.282222) q[0];
rz(1.1497633) q[1];
sx q[1];
rz(-1.7825922) q[1];
sx q[1];
rz(-1.5001635) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.951913) q[0];
sx q[0];
rz(-1.743058) q[0];
sx q[0];
rz(0.36393117) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9740536) q[2];
sx q[2];
rz(-2.5563498) q[2];
sx q[2];
rz(-0.61589059) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0974429) q[1];
sx q[1];
rz(-1.2359834) q[1];
sx q[1];
rz(-2.8756932) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30434609) q[3];
sx q[3];
rz(-1.658421) q[3];
sx q[3];
rz(0.1803785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.942261) q[2];
sx q[2];
rz(-1.4456238) q[2];
sx q[2];
rz(-1.65421) q[2];
rz(-2.2516294) q[3];
sx q[3];
rz(-1.7421937) q[3];
sx q[3];
rz(-0.14931211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99059659) q[0];
sx q[0];
rz(-0.52495933) q[0];
sx q[0];
rz(1.6816444) q[0];
rz(-1.2851985) q[1];
sx q[1];
rz(-1.3672914) q[1];
sx q[1];
rz(-0.45305124) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.995034) q[0];
sx q[0];
rz(-1.5164598) q[0];
sx q[0];
rz(1.8749692) q[0];
x q[1];
rz(1.4905246) q[2];
sx q[2];
rz(-2.5577214) q[2];
sx q[2];
rz(-0.50498) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1374319) q[1];
sx q[1];
rz(-0.42037005) q[1];
sx q[1];
rz(2.67291) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3052081) q[3];
sx q[3];
rz(-2.5646113) q[3];
sx q[3];
rz(-1.0671237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.65752658) q[2];
sx q[2];
rz(-2.6675197) q[2];
sx q[2];
rz(1.5566114) q[2];
rz(-1.5629684) q[3];
sx q[3];
rz(-1.862674) q[3];
sx q[3];
rz(-1.0274308) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9180561) q[0];
sx q[0];
rz(-2.1503088) q[0];
sx q[0];
rz(1.9687442) q[0];
rz(-0.46074834) q[1];
sx q[1];
rz(-1.8811503) q[1];
sx q[1];
rz(1.4985098) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39678412) q[0];
sx q[0];
rz(-0.61273324) q[0];
sx q[0];
rz(1.5721129) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4732359) q[2];
sx q[2];
rz(-0.41497013) q[2];
sx q[2];
rz(2.9573832) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.071049404) q[1];
sx q[1];
rz(-0.26916801) q[1];
sx q[1];
rz(1.7885923) q[1];
rz(-pi) q[2];
rz(-3.0027578) q[3];
sx q[3];
rz(-2.2679726) q[3];
sx q[3];
rz(2.4133854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1341165) q[2];
sx q[2];
rz(-1.4223998) q[2];
sx q[2];
rz(0.22809347) q[2];
rz(-0.54276931) q[3];
sx q[3];
rz(-1.0017064) q[3];
sx q[3];
rz(-2.2293279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-3.0246564) q[0];
sx q[0];
rz(-1.0448562) q[0];
sx q[0];
rz(0.60633099) q[0];
rz(1.2888651) q[1];
sx q[1];
rz(-1.6090798) q[1];
sx q[1];
rz(-0.85561633) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80911773) q[0];
sx q[0];
rz(-2.1590589) q[0];
sx q[0];
rz(2.0659129) q[0];
rz(-1.6560203) q[2];
sx q[2];
rz(-2.3422675) q[2];
sx q[2];
rz(2.8003789) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2825477) q[1];
sx q[1];
rz(-0.63524109) q[1];
sx q[1];
rz(-0.35104351) q[1];
x q[2];
rz(-1.0904013) q[3];
sx q[3];
rz(-1.3452936) q[3];
sx q[3];
rz(0.71798872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3768846) q[2];
sx q[2];
rz(-0.43262425) q[2];
sx q[2];
rz(-3.063859) q[2];
rz(-3.0149095) q[3];
sx q[3];
rz(-2.099791) q[3];
sx q[3];
rz(-2.3104987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7938101) q[0];
sx q[0];
rz(-0.3599444) q[0];
sx q[0];
rz(-2.1600294) q[0];
rz(1.7836102) q[1];
sx q[1];
rz(-0.49109083) q[1];
sx q[1];
rz(0.29409274) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2002073) q[0];
sx q[0];
rz(-1.6305171) q[0];
sx q[0];
rz(-2.8325547) q[0];
rz(2.5035759) q[2];
sx q[2];
rz(-2.6891481) q[2];
sx q[2];
rz(-2.5663301) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5433054) q[1];
sx q[1];
rz(-1.2906055) q[1];
sx q[1];
rz(2.2375475) q[1];
x q[2];
rz(-0.050042583) q[3];
sx q[3];
rz(-1.4980157) q[3];
sx q[3];
rz(0.67228729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0391417) q[2];
sx q[2];
rz(-0.65639085) q[2];
sx q[2];
rz(1.0605109) q[2];
rz(2.5833526) q[3];
sx q[3];
rz(-0.57359901) q[3];
sx q[3];
rz(-0.019088117) q[3];
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
rz(-1.7790262) q[0];
sx q[0];
rz(-2.9988852) q[0];
sx q[0];
rz(-2.3574164) q[0];
rz(-1.342429) q[1];
sx q[1];
rz(-1.8524086) q[1];
sx q[1];
rz(2.4450891) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9331037) q[0];
sx q[0];
rz(-1.3580926) q[0];
sx q[0];
rz(-1.3377633) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4694211) q[2];
sx q[2];
rz(-1.1706588) q[2];
sx q[2];
rz(-0.53022417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3236774) q[1];
sx q[1];
rz(-1.9888048) q[1];
sx q[1];
rz(0.54115414) q[1];
rz(1.7197929) q[3];
sx q[3];
rz(-1.7033352) q[3];
sx q[3];
rz(0.46035351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33783087) q[2];
sx q[2];
rz(-0.75118128) q[2];
sx q[2];
rz(-1.2639698) q[2];
rz(1.3761282) q[3];
sx q[3];
rz(-0.91457808) q[3];
sx q[3];
rz(-1.5553364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92397583) q[0];
sx q[0];
rz(-0.042984977) q[0];
sx q[0];
rz(1.6643583) q[0];
rz(0.90786511) q[1];
sx q[1];
rz(-1.6198747) q[1];
sx q[1];
rz(2.1715865) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15031397) q[0];
sx q[0];
rz(-1.3320001) q[0];
sx q[0];
rz(2.1697609) q[0];
rz(-pi) q[1];
rz(1.1437505) q[2];
sx q[2];
rz(-2.5055024) q[2];
sx q[2];
rz(1.191312) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.86783389) q[1];
sx q[1];
rz(-0.9504488) q[1];
sx q[1];
rz(-1.7628487) q[1];
rz(-pi) q[2];
rz(-2.41955) q[3];
sx q[3];
rz(-2.0835428) q[3];
sx q[3];
rz(1.3861314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9475391) q[2];
sx q[2];
rz(-0.87475646) q[2];
sx q[2];
rz(0.27210316) q[2];
rz(-0.84723204) q[3];
sx q[3];
rz(-2.6341485) q[3];
sx q[3];
rz(2.0525172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.465441) q[0];
sx q[0];
rz(-1.1347102) q[0];
sx q[0];
rz(1.4665428) q[0];
rz(2.8027986) q[1];
sx q[1];
rz(-1.6962961) q[1];
sx q[1];
rz(1.2512339) q[1];
rz(2.4215578) q[2];
sx q[2];
rz(-1.3007813) q[2];
sx q[2];
rz(-0.5310975) q[2];
rz(0.083214464) q[3];
sx q[3];
rz(-1.1490001) q[3];
sx q[3];
rz(0.84925539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
