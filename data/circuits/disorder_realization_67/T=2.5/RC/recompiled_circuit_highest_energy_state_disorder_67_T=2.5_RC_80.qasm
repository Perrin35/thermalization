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
rz(-0.026365658) q[0];
rz(0.68323505) q[1];
sx q[1];
rz(4.9024138) q[1];
sx q[1];
rz(9.496357) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3535093) q[0];
sx q[0];
rz(-3.1147163) q[0];
sx q[0];
rz(0.36969097) q[0];
rz(1.8060922) q[2];
sx q[2];
rz(-1.9517731) q[2];
sx q[2];
rz(1.345696) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72373448) q[1];
sx q[1];
rz(-1.5520387) q[1];
sx q[1];
rz(0.0049519227) q[1];
x q[2];
rz(-0.79238331) q[3];
sx q[3];
rz(-0.43607831) q[3];
sx q[3];
rz(0.60286264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8611531) q[2];
sx q[2];
rz(-2.4344378) q[2];
sx q[2];
rz(-0.46671483) q[2];
rz(2.6672582) q[3];
sx q[3];
rz(-3.1204087) q[3];
sx q[3];
rz(-3.1202313) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83653432) q[0];
sx q[0];
rz(-0.49764043) q[0];
sx q[0];
rz(-3.1217788) q[0];
rz(-1.548798) q[1];
sx q[1];
rz(-2.9217547) q[1];
sx q[1];
rz(-1.4733431) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9813741) q[0];
sx q[0];
rz(-0.88079534) q[0];
sx q[0];
rz(-0.86440683) q[0];
rz(-pi) q[1];
rz(-2.8823518) q[2];
sx q[2];
rz(-1.8564285) q[2];
sx q[2];
rz(2.223857) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87621385) q[1];
sx q[1];
rz(-2.9326831) q[1];
sx q[1];
rz(1.9494328) q[1];
rz(-pi) q[2];
rz(-0.83062828) q[3];
sx q[3];
rz(-2.1601956) q[3];
sx q[3];
rz(2.7089861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3026498) q[2];
sx q[2];
rz(-2.5431716) q[2];
sx q[2];
rz(-1.8518651) q[2];
rz(1.2368115) q[3];
sx q[3];
rz(-0.33331063) q[3];
sx q[3];
rz(-0.70241565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.1692093) q[0];
sx q[0];
rz(-1.9820259) q[0];
sx q[0];
rz(-1.5706536) q[0];
rz(-1.6698569) q[1];
sx q[1];
rz(-1.6153299) q[1];
sx q[1];
rz(0.45438802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1412068) q[0];
sx q[0];
rz(-2.1659746) q[0];
sx q[0];
rz(-2.9808874) q[0];
rz(-pi) q[1];
rz(-1.4287895) q[2];
sx q[2];
rz(-1.5446086) q[2];
sx q[2];
rz(-2.533503) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4799166) q[1];
sx q[1];
rz(-3.0212086) q[1];
sx q[1];
rz(-1.7604339) q[1];
rz(-2.2265395) q[3];
sx q[3];
rz(-2.1106535) q[3];
sx q[3];
rz(-1.2815544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4001974) q[2];
sx q[2];
rz(-0.016308451) q[2];
sx q[2];
rz(2.7941217) q[2];
rz(-2.6853284) q[3];
sx q[3];
rz(-0.01472344) q[3];
sx q[3];
rz(2.1485476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149813) q[0];
sx q[0];
rz(-1.9009637) q[0];
sx q[0];
rz(-1.4062784) q[0];
rz(-2.6982488) q[1];
sx q[1];
rz(-2.1163546) q[1];
sx q[1];
rz(-1.5700856) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3298108) q[0];
sx q[0];
rz(-0.64231163) q[0];
sx q[0];
rz(-0.091600939) q[0];
rz(-pi) q[1];
x q[1];
rz(0.069025682) q[2];
sx q[2];
rz(-1.6352425) q[2];
sx q[2];
rz(0.2914337) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52165495) q[1];
sx q[1];
rz(-1.8490377) q[1];
sx q[1];
rz(3.131429) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27914234) q[3];
sx q[3];
rz(-1.9873535) q[3];
sx q[3];
rz(-3.109351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7683679) q[2];
sx q[2];
rz(-0.39373213) q[2];
sx q[2];
rz(3.0623398) q[2];
rz(1.9521889) q[3];
sx q[3];
rz(-1.7704084) q[3];
sx q[3];
rz(1.5090548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4382512) q[0];
sx q[0];
rz(-0.51333135) q[0];
sx q[0];
rz(-0.80605036) q[0];
rz(2.3015859) q[1];
sx q[1];
rz(-0.012931074) q[1];
sx q[1];
rz(0.77617019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38473338) q[0];
sx q[0];
rz(-0.98664588) q[0];
sx q[0];
rz(-1.4136876) q[0];
rz(-pi) q[1];
rz(-1.569328) q[2];
sx q[2];
rz(-1.5812354) q[2];
sx q[2];
rz(-0.28077747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.35368109) q[1];
sx q[1];
rz(-3.0056142) q[1];
sx q[1];
rz(-0.087894364) q[1];
rz(-pi) q[2];
rz(-2.8819894) q[3];
sx q[3];
rz(-1.6791376) q[3];
sx q[3];
rz(-1.3110127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.025803056) q[2];
sx q[2];
rz(-1.5805406) q[2];
sx q[2];
rz(0.72581327) q[2];
rz(-0.18802655) q[3];
sx q[3];
rz(-0.060421061) q[3];
sx q[3];
rz(2.4316725) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96227729) q[0];
sx q[0];
rz(-2.582452) q[0];
sx q[0];
rz(0.54139262) q[0];
rz(2.9560282) q[1];
sx q[1];
rz(-1.5488397) q[1];
sx q[1];
rz(0.12841368) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5914766) q[0];
sx q[0];
rz(-1.8719561) q[0];
sx q[0];
rz(0.37806435) q[0];
rz(1.6921243) q[2];
sx q[2];
rz(-1.5660966) q[2];
sx q[2];
rz(0.0041065816) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7306108) q[1];
sx q[1];
rz(-1.1514947) q[1];
sx q[1];
rz(1.4813408) q[1];
rz(-pi) q[2];
rz(-1.7592247) q[3];
sx q[3];
rz(-1.5811833) q[3];
sx q[3];
rz(0.56558181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.386117) q[2];
sx q[2];
rz(-3.0839034) q[2];
sx q[2];
rz(2.3003787) q[2];
rz(2.9403213) q[3];
sx q[3];
rz(-1.6095251) q[3];
sx q[3];
rz(2.8819528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5801308) q[0];
sx q[0];
rz(-2.353297) q[0];
sx q[0];
rz(-1.5808251) q[0];
rz(-0.67281094) q[1];
sx q[1];
rz(-1.4038059) q[1];
sx q[1];
rz(0.028884551) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2530841) q[0];
sx q[0];
rz(-0.38034236) q[0];
sx q[0];
rz(-2.4230291) q[0];
x q[1];
rz(-2.6631764) q[2];
sx q[2];
rz(-1.3923651) q[2];
sx q[2];
rz(-2.894424) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1950732) q[1];
sx q[1];
rz(-1.9819489) q[1];
sx q[1];
rz(-2.8821936) q[1];
rz(-pi) q[2];
rz(-2.1013124) q[3];
sx q[3];
rz(-0.86581745) q[3];
sx q[3];
rz(1.353273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9081356) q[2];
sx q[2];
rz(-0.59671777) q[2];
sx q[2];
rz(-0.92823589) q[2];
rz(-0.2975896) q[3];
sx q[3];
rz(-2.987515) q[3];
sx q[3];
rz(-2.2969864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069227844) q[0];
sx q[0];
rz(-2.897825) q[0];
sx q[0];
rz(0.099076554) q[0];
rz(0.98723269) q[1];
sx q[1];
rz(-1.8376553) q[1];
sx q[1];
rz(-2.6027021) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9469556) q[0];
sx q[0];
rz(-1.4879003) q[0];
sx q[0];
rz(-3.1350053) q[0];
x q[1];
rz(1.107223) q[2];
sx q[2];
rz(-3.0863783) q[2];
sx q[2];
rz(-0.56978536) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3618631) q[1];
sx q[1];
rz(-1.9302082) q[1];
sx q[1];
rz(0.56198175) q[1];
x q[2];
rz(1.6383) q[3];
sx q[3];
rz(-1.5482404) q[3];
sx q[3];
rz(-1.8543138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99523669) q[2];
sx q[2];
rz(-0.015492798) q[2];
sx q[2];
rz(-2.8323979) q[2];
rz(2.501798) q[3];
sx q[3];
rz(-0.00034172405) q[3];
sx q[3];
rz(3.0777847) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30534202) q[0];
sx q[0];
rz(-2.5449365) q[0];
sx q[0];
rz(-0.054656595) q[0];
rz(-1.1326185) q[1];
sx q[1];
rz(-1.9839958) q[1];
sx q[1];
rz(-1.2861015) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085145935) q[0];
sx q[0];
rz(-2.9582199) q[0];
sx q[0];
rz(1.3297179) q[0];
rz(3.1007721) q[2];
sx q[2];
rz(-1.5287182) q[2];
sx q[2];
rz(0.034860858) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1117592) q[1];
sx q[1];
rz(-2.8984424) q[1];
sx q[1];
rz(-2.7584226) q[1];
rz(2.0580106) q[3];
sx q[3];
rz(-1.6346674) q[3];
sx q[3];
rz(2.4364249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9201811) q[2];
sx q[2];
rz(-2.5765918) q[2];
sx q[2];
rz(2.0378225) q[2];
rz(-1.6086027) q[3];
sx q[3];
rz(-3.0975603) q[3];
sx q[3];
rz(0.65582961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.00103818) q[0];
sx q[0];
rz(-0.1796722) q[0];
sx q[0];
rz(3.1371064) q[0];
rz(1.5525612) q[1];
sx q[1];
rz(-1.4493425) q[1];
sx q[1];
rz(-3.0844614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.321314) q[0];
sx q[0];
rz(-1.4749267) q[0];
sx q[0];
rz(1.7509599) q[0];
rz(-pi) q[1];
rz(-0.22259076) q[2];
sx q[2];
rz(-1.711297) q[2];
sx q[2];
rz(-2.8369571) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.58913402) q[1];
sx q[1];
rz(-1.4499364) q[1];
sx q[1];
rz(-2.2446988) q[1];
x q[2];
rz(-0.1520429) q[3];
sx q[3];
rz(-0.83629823) q[3];
sx q[3];
rz(1.9805857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39157465) q[2];
sx q[2];
rz(-0.027316814) q[2];
sx q[2];
rz(2.2726783) q[2];
rz(-2.1598375) q[3];
sx q[3];
rz(-0.029564094) q[3];
sx q[3];
rz(0.50299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.045573087) q[0];
sx q[0];
rz(-1.4843142) q[0];
sx q[0];
rz(1.6577161) q[0];
rz(-2.6847196) q[1];
sx q[1];
rz(-2.9869106) q[1];
sx q[1];
rz(3.0966495) q[1];
rz(1.7560319) q[2];
sx q[2];
rz(-2.3262263) q[2];
sx q[2];
rz(-0.066166886) q[2];
rz(-1.4275238) q[3];
sx q[3];
rz(-2.2289056) q[3];
sx q[3];
rz(-3.0994305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
