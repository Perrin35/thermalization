OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93513918) q[0];
sx q[0];
rz(3.9225188) q[0];
sx q[0];
rz(9.6315686) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(0.18016711) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8040112) q[0];
sx q[0];
rz(-1.1085701) q[0];
sx q[0];
rz(1.1929212) q[0];
rz(-1.4397058) q[2];
sx q[2];
rz(-1.0616454) q[2];
sx q[2];
rz(0.57123643) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.78046103) q[1];
sx q[1];
rz(-1.4968922) q[1];
sx q[1];
rz(-0.88688811) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7546685) q[3];
sx q[3];
rz(-2.2443218) q[3];
sx q[3];
rz(-0.20007381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41123286) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(-1.8475378) q[2];
rz(0.40575746) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(0.40677235) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2186573) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(-0.12582114) q[0];
rz(-2.3361092) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(-1.3719826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3208647) q[0];
sx q[0];
rz(-0.67824368) q[0];
sx q[0];
rz(0.098537785) q[0];
x q[1];
rz(-1.6945018) q[2];
sx q[2];
rz(-2.1134085) q[2];
sx q[2];
rz(-1.5109085) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6276715) q[1];
sx q[1];
rz(-1.5253592) q[1];
sx q[1];
rz(-2.8915878) q[1];
x q[2];
rz(0.45687859) q[3];
sx q[3];
rz(-2.6279454) q[3];
sx q[3];
rz(1.160624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25386086) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(2.1453693) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(-2.2912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4662194) q[0];
sx q[0];
rz(-0.96005625) q[0];
sx q[0];
rz(-1.0590142) q[0];
rz(-1.9937218) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(-1.0645197) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36944593) q[0];
sx q[0];
rz(-1.1755953) q[0];
sx q[0];
rz(-0.69338436) q[0];
rz(-1.9532922) q[2];
sx q[2];
rz(-1.7171535) q[2];
sx q[2];
rz(2.086703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3036084) q[1];
sx q[1];
rz(-1.8808937) q[1];
sx q[1];
rz(3.1293948) q[1];
x q[2];
rz(-2.2952609) q[3];
sx q[3];
rz(-1.4457821) q[3];
sx q[3];
rz(2.4785329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.069783) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(-2.8783669) q[2];
rz(-2.0227382) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(-0.81937218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333106) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(-1.7472965) q[0];
rz(1.0385723) q[1];
sx q[1];
rz(-1.4515406) q[1];
sx q[1];
rz(-1.0669473) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79558668) q[0];
sx q[0];
rz(-0.53115244) q[0];
sx q[0];
rz(1.4941494) q[0];
rz(-pi) q[1];
rz(2.8720886) q[2];
sx q[2];
rz(-1.6490893) q[2];
sx q[2];
rz(1.6646202) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9001138) q[1];
sx q[1];
rz(-2.334324) q[1];
sx q[1];
rz(0.57358731) q[1];
rz(-pi) q[2];
rz(-2.0526485) q[3];
sx q[3];
rz(-1.7715766) q[3];
sx q[3];
rz(2.5364385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84714326) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(-0.28953826) q[2];
rz(0.52044049) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(0.7590487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(-0.83475137) q[0];
rz(-1.2443776) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(1.429819) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0833417) q[0];
sx q[0];
rz(-1.6334851) q[0];
sx q[0];
rz(-1.5682674) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3313815) q[2];
sx q[2];
rz(-1.8481701) q[2];
sx q[2];
rz(-1.1053567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.012638481) q[1];
sx q[1];
rz(-1.5938938) q[1];
sx q[1];
rz(0.79146339) q[1];
rz(-pi) q[2];
rz(1.7992713) q[3];
sx q[3];
rz(-0.67111525) q[3];
sx q[3];
rz(1.6959147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66118583) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(-1.2832114) q[2];
rz(-3.0094106) q[3];
sx q[3];
rz(-2.81288) q[3];
sx q[3];
rz(-2.8018518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(-0.47750372) q[0];
rz(1.6409138) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(-2.5710411) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084224) q[0];
sx q[0];
rz(-2.6714984) q[0];
sx q[0];
rz(0.8870468) q[0];
rz(-2.4628377) q[2];
sx q[2];
rz(-2.5828913) q[2];
sx q[2];
rz(2.9847381) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8385982) q[1];
sx q[1];
rz(-1.0698776) q[1];
sx q[1];
rz(-2.3272728) q[1];
rz(-pi) q[2];
rz(1.5485351) q[3];
sx q[3];
rz(-2.2189757) q[3];
sx q[3];
rz(2.9472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.70696124) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(0.79664191) q[2];
rz(-0.32026511) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(-1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758078) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(2.7835223) q[0];
rz(-2.8529196) q[1];
sx q[1];
rz(-0.52983785) q[1];
sx q[1];
rz(-1.3105062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.34328) q[0];
sx q[0];
rz(-0.57050059) q[0];
sx q[0];
rz(-2.7891854) q[0];
x q[1];
rz(-2.481776) q[2];
sx q[2];
rz(-0.81507999) q[2];
sx q[2];
rz(1.3224524) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6984182) q[1];
sx q[1];
rz(-2.0768754) q[1];
sx q[1];
rz(-2.0884104) q[1];
rz(-pi) q[2];
rz(0.88498022) q[3];
sx q[3];
rz(-1.0776057) q[3];
sx q[3];
rz(-2.5964338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8075809) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(0.99747783) q[2];
rz(-0.57972646) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8758133) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(0.34061256) q[0];
rz(-1.9050725) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(-1.9514726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1163568) q[0];
sx q[0];
rz(-1.802117) q[0];
sx q[0];
rz(0.063409253) q[0];
rz(-0.868452) q[2];
sx q[2];
rz(-3.0358899) q[2];
sx q[2];
rz(0.58123523) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78672532) q[1];
sx q[1];
rz(-2.7495972) q[1];
sx q[1];
rz(0.74704945) q[1];
rz(-0.45300071) q[3];
sx q[3];
rz(-2.0815597) q[3];
sx q[3];
rz(2.9651027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3999346) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(2.7271872) q[2];
rz(1.7587781) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8906616) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(-2.9898306) q[0];
rz(-1.4258619) q[1];
sx q[1];
rz(-1.3769923) q[1];
sx q[1];
rz(0.94917667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28578368) q[0];
sx q[0];
rz(-1.7109509) q[0];
sx q[0];
rz(2.6967718) q[0];
rz(-pi) q[1];
rz(-0.62990909) q[2];
sx q[2];
rz(-1.7498778) q[2];
sx q[2];
rz(-1.5392898) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9526082) q[1];
sx q[1];
rz(-2.4338255) q[1];
sx q[1];
rz(-2.7680552) q[1];
x q[2];
rz(-0.98032326) q[3];
sx q[3];
rz(-1.1579517) q[3];
sx q[3];
rz(2.1291898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40733797) q[2];
sx q[2];
rz(-0.39525017) q[2];
sx q[2];
rz(-2.7091743) q[2];
rz(1.6010823) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0712873) q[0];
sx q[0];
rz(-0.27305958) q[0];
sx q[0];
rz(-0.37316698) q[0];
rz(-0.35692731) q[1];
sx q[1];
rz(-0.86078763) q[1];
sx q[1];
rz(-2.4180791) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31691566) q[0];
sx q[0];
rz(-0.74130374) q[0];
sx q[0];
rz(2.9497428) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4805484) q[2];
sx q[2];
rz(-0.89333488) q[2];
sx q[2];
rz(-1.2573164) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8404624) q[1];
sx q[1];
rz(-2.3311619) q[1];
sx q[1];
rz(0.37709548) q[1];
x q[2];
rz(-0.27042737) q[3];
sx q[3];
rz(-1.575483) q[3];
sx q[3];
rz(3.1386496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2202806) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(0.55345654) q[2];
rz(-0.71183318) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.9217459) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(-0.28221054) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(-1.1214592) q[2];
sx q[2];
rz(-1.1171787) q[2];
sx q[2];
rz(1.6906307) q[2];
rz(2.5440352) q[3];
sx q[3];
rz(-0.47511027) q[3];
sx q[3];
rz(-2.6078754) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
