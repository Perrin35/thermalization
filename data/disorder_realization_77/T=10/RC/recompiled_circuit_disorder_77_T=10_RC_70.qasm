OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(-2.9736829) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(-2.8462703) q[1];
sx q[1];
rz(0.056161031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18626285) q[0];
sx q[0];
rz(-1.7832527) q[0];
sx q[0];
rz(-2.0583378) q[0];
rz(-pi) q[1];
rz(-1.2916318) q[2];
sx q[2];
rz(-2.4764875) q[2];
sx q[2];
rz(-1.4814188) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.72370428) q[1];
sx q[1];
rz(-1.0350409) q[1];
sx q[1];
rz(-2.8292311) q[1];
x q[2];
rz(-1.6151186) q[3];
sx q[3];
rz(-1.830415) q[3];
sx q[3];
rz(-2.0093902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(2.7089233) q[2];
rz(-1.9487322) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(-1.8288076) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(1.9899433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5702471) q[0];
sx q[0];
rz(-0.76128188) q[0];
sx q[0];
rz(-1.135457) q[0];
rz(0.67302455) q[2];
sx q[2];
rz(-2.4403799) q[2];
sx q[2];
rz(2.6075624) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.093106) q[1];
sx q[1];
rz(-1.9890607) q[1];
sx q[1];
rz(-2.6622245) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1851951) q[3];
sx q[3];
rz(-2.0413627) q[3];
sx q[3];
rz(3.1009931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(0.22182626) q[2];
rz(-2.7644073) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8310228) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(3.1047399) q[0];
rz(-0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(3.085014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.409257) q[0];
sx q[0];
rz(-1.1261254) q[0];
sx q[0];
rz(2.1952941) q[0];
rz(-1.4372196) q[2];
sx q[2];
rz(-1.3165511) q[2];
sx q[2];
rz(-0.1453407) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50443447) q[1];
sx q[1];
rz(-1.8706026) q[1];
sx q[1];
rz(-1.7791041) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64485456) q[3];
sx q[3];
rz(-1.8861024) q[3];
sx q[3];
rz(-0.23526084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-2.2154714) q[2];
rz(-0.55666322) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(-1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2264003) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(2.6779968) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6620561) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(2.7583073) q[0];
x q[1];
rz(0.17642994) q[2];
sx q[2];
rz(-2.8039805) q[2];
sx q[2];
rz(-2.7531429) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78471781) q[1];
sx q[1];
rz(-2.3349635) q[1];
sx q[1];
rz(2.9033317) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0158402) q[3];
sx q[3];
rz(-1.1557126) q[3];
sx q[3];
rz(-2.2500028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7745557) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(-3.0920933) q[2];
rz(-0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(0.4822576) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5807242) q[0];
sx q[0];
rz(-1.3618016) q[0];
sx q[0];
rz(-3.0955549) q[0];
x q[1];
rz(-2.7371251) q[2];
sx q[2];
rz(-2.3239115) q[2];
sx q[2];
rz(2.5596465) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3790834) q[1];
sx q[1];
rz(-0.06113872) q[1];
sx q[1];
rz(1.5361384) q[1];
rz(-pi) q[2];
rz(-1.3145507) q[3];
sx q[3];
rz(-1.991193) q[3];
sx q[3];
rz(2.8375569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7047983) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(2.5571402) q[0];
rz(2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(-0.054919682) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7391101) q[0];
sx q[0];
rz(-1.4882898) q[0];
sx q[0];
rz(1.8399747) q[0];
x q[1];
rz(1.6266277) q[2];
sx q[2];
rz(-2.8179114) q[2];
sx q[2];
rz(-1.2602381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1153032) q[1];
sx q[1];
rz(-2.4541828) q[1];
sx q[1];
rz(0.61647146) q[1];
x q[2];
rz(-2.1543598) q[3];
sx q[3];
rz(-1.9658486) q[3];
sx q[3];
rz(2.6810255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(-1.0167271) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(-1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-2.8570535) q[0];
rz(0.94447213) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(-0.91032666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9757662) q[0];
sx q[0];
rz(-1.5162139) q[0];
sx q[0];
rz(-1.635701) q[0];
rz(-0.89332135) q[2];
sx q[2];
rz(-2.6258694) q[2];
sx q[2];
rz(2.233778) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7200304) q[1];
sx q[1];
rz(-1.7517462) q[1];
sx q[1];
rz(-0.49444316) q[1];
rz(-pi) q[2];
rz(-0.72083731) q[3];
sx q[3];
rz(-1.1158873) q[3];
sx q[3];
rz(2.7618559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(2.2195623) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(-2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2475125) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-0.12938736) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(-0.30050373) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43113118) q[0];
sx q[0];
rz(-2.3346402) q[0];
sx q[0];
rz(-2.0956844) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6408765) q[2];
sx q[2];
rz(-2.2705728) q[2];
sx q[2];
rz(-1.7388294) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1860736) q[1];
sx q[1];
rz(-1.3961853) q[1];
sx q[1];
rz(-2.7956635) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81340202) q[3];
sx q[3];
rz(-1.3374995) q[3];
sx q[3];
rz(-0.63026159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5552716) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(0.55109763) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-0.12776275) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(0.75884563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7153873) q[0];
sx q[0];
rz(-1.1466768) q[0];
sx q[0];
rz(3.12294) q[0];
rz(-pi) q[1];
rz(-0.60894572) q[2];
sx q[2];
rz(-1.4119838) q[2];
sx q[2];
rz(1.9132523) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85941852) q[1];
sx q[1];
rz(-2.4255883) q[1];
sx q[1];
rz(2.9031258) q[1];
rz(-3.127029) q[3];
sx q[3];
rz(-0.90523883) q[3];
sx q[3];
rz(-1.927782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1252497) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(0.49003595) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-0.075335659) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(-2.5316701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1760575) q[0];
sx q[0];
rz(-2.3241204) q[0];
sx q[0];
rz(-2.7479991) q[0];
rz(2.7650325) q[2];
sx q[2];
rz(-1.8278404) q[2];
sx q[2];
rz(3.0388447) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.016644195) q[1];
sx q[1];
rz(-0.84608191) q[1];
sx q[1];
rz(-1.8925341) q[1];
x q[2];
rz(0.742357) q[3];
sx q[3];
rz(-1.9543813) q[3];
sx q[3];
rz(-2.50768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23218368) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(0.71371901) q[2];
rz(-0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(-2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80355766) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(-3.1109839) q[2];
sx q[2];
rz(-1.7666713) q[2];
sx q[2];
rz(-0.91790744) q[2];
rz(-1.0999023) q[3];
sx q[3];
rz(-1.7970016) q[3];
sx q[3];
rz(-2.5947528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
