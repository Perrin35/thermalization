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
rz(3.115227) q[0];
rz(0.68323505) q[1];
sx q[1];
rz(-1.3807715) q[1];
sx q[1];
rz(0.071579054) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.723322) q[0];
sx q[0];
rz(-1.5958565) q[0];
sx q[0];
rz(-1.5805095) q[0];
rz(1.3355005) q[2];
sx q[2];
rz(-1.9517731) q[2];
sx q[2];
rz(-1.345696) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4178582) q[1];
sx q[1];
rz(-1.5520387) q[1];
sx q[1];
rz(3.1366407) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31622131) q[3];
sx q[3];
rz(-1.8762778) q[3];
sx q[3];
rz(-0.22462928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2804395) q[2];
sx q[2];
rz(-2.4344378) q[2];
sx q[2];
rz(-2.6748778) q[2];
rz(2.6672582) q[3];
sx q[3];
rz(-3.1204087) q[3];
sx q[3];
rz(-3.1202313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83653432) q[0];
sx q[0];
rz(-0.49764043) q[0];
sx q[0];
rz(0.019813892) q[0];
rz(1.5927947) q[1];
sx q[1];
rz(-2.9217547) q[1];
sx q[1];
rz(-1.4733431) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9813741) q[0];
sx q[0];
rz(-2.2607973) q[0];
sx q[0];
rz(-2.2771858) q[0];
rz(-pi) q[1];
rz(0.25924087) q[2];
sx q[2];
rz(-1.8564285) q[2];
sx q[2];
rz(2.223857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.87621385) q[1];
sx q[1];
rz(-0.20890954) q[1];
sx q[1];
rz(-1.9494328) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4056689) q[3];
sx q[3];
rz(-0.97565996) q[3];
sx q[3];
rz(1.6079966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8389429) q[2];
sx q[2];
rz(-0.5984211) q[2];
sx q[2];
rz(-1.8518651) q[2];
rz(-1.9047811) q[3];
sx q[3];
rz(-0.33331063) q[3];
sx q[3];
rz(2.439177) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1692093) q[0];
sx q[0];
rz(-1.1595668) q[0];
sx q[0];
rz(1.5709391) q[0];
rz(-1.4717357) q[1];
sx q[1];
rz(-1.5262628) q[1];
sx q[1];
rz(-2.6872046) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1412068) q[0];
sx q[0];
rz(-2.1659746) q[0];
sx q[0];
rz(-2.9808874) q[0];
rz(-pi) q[1];
rz(1.4287895) q[2];
sx q[2];
rz(-1.5446086) q[2];
sx q[2];
rz(-0.60808966) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4799166) q[1];
sx q[1];
rz(-0.12038409) q[1];
sx q[1];
rz(-1.7604339) q[1];
rz(-pi) q[2];
rz(-2.4942355) q[3];
sx q[3];
rz(-1.0203927) q[3];
sx q[3];
rz(0.087322012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74139524) q[2];
sx q[2];
rz(-0.016308451) q[2];
sx q[2];
rz(2.7941217) q[2];
rz(2.6853284) q[3];
sx q[3];
rz(-3.1268692) q[3];
sx q[3];
rz(-0.99304503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149813) q[0];
sx q[0];
rz(-1.9009637) q[0];
sx q[0];
rz(1.7353143) q[0];
rz(-2.6982488) q[1];
sx q[1];
rz(-1.025238) q[1];
sx q[1];
rz(1.5700856) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8271874) q[0];
sx q[0];
rz(-1.6256204) q[0];
sx q[0];
rz(2.501295) q[0];
rz(-0.75218828) q[2];
sx q[2];
rz(-3.0471932) q[2];
sx q[2];
rz(2.029325) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0896596) q[1];
sx q[1];
rz(-1.580569) q[1];
sx q[1];
rz(1.2925413) q[1];
x q[2];
rz(-1.1394386) q[3];
sx q[3];
rz(-1.8255116) q[3];
sx q[3];
rz(1.4875808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37322474) q[2];
sx q[2];
rz(-0.39373213) q[2];
sx q[2];
rz(-3.0623398) q[2];
rz(1.9521889) q[3];
sx q[3];
rz(-1.7704084) q[3];
sx q[3];
rz(-1.6325379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-3.1286616) q[1];
sx q[1];
rz(2.3654225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66444976) q[0];
sx q[0];
rz(-0.60252554) q[0];
sx q[0];
rz(0.23238927) q[0];
x q[1];
rz(3.0018535) q[2];
sx q[2];
rz(-0.010541803) q[2];
sx q[2];
rz(0.14103061) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8373903) q[1];
sx q[1];
rz(-1.5826962) q[1];
sx q[1];
rz(0.13545998) q[1];
rz(-pi) q[2];
rz(1.4587298) q[3];
sx q[3];
rz(-1.8288426) q[3];
sx q[3];
rz(0.28849453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1157896) q[2];
sx q[2];
rz(-1.5610521) q[2];
sx q[2];
rz(0.72581327) q[2];
rz(-0.18802655) q[3];
sx q[3];
rz(-3.0811716) q[3];
sx q[3];
rz(-2.4316725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1793154) q[0];
sx q[0];
rz(-0.55914068) q[0];
sx q[0];
rz(2.6002) q[0];
rz(0.18556449) q[1];
sx q[1];
rz(-1.5927529) q[1];
sx q[1];
rz(-3.013179) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62080431) q[0];
sx q[0];
rz(-2.6628011) q[0];
sx q[0];
rz(-0.6995246) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0047345078) q[2];
sx q[2];
rz(-1.4494697) q[2];
sx q[2];
rz(-1.5743299) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12331387) q[1];
sx q[1];
rz(-1.4891081) q[1];
sx q[1];
rz(-2.7207992) q[1];
rz(-pi) q[2];
rz(-0.010574118) q[3];
sx q[3];
rz(-1.3823783) q[3];
sx q[3];
rz(2.1343975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7554756) q[2];
sx q[2];
rz(-3.0839034) q[2];
sx q[2];
rz(0.84121394) q[2];
rz(2.9403213) q[3];
sx q[3];
rz(-1.6095251) q[3];
sx q[3];
rz(-0.25963983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5614618) q[0];
sx q[0];
rz(-0.78829563) q[0];
sx q[0];
rz(-1.5607675) q[0];
rz(-0.67281094) q[1];
sx q[1];
rz(-1.4038059) q[1];
sx q[1];
rz(-3.1127081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7772693) q[0];
sx q[0];
rz(-1.3239081) q[0];
sx q[0];
rz(2.8492575) q[0];
x q[1];
rz(-2.7682226) q[2];
sx q[2];
rz(-0.50818102) q[2];
sx q[2];
rz(2.1477107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.72994443) q[1];
sx q[1];
rz(-1.333451) q[1];
sx q[1];
rz(-1.1470331) q[1];
x q[2];
rz(-2.6050909) q[3];
sx q[3];
rz(-2.287545) q[3];
sx q[3];
rz(-2.0889747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2334571) q[2];
sx q[2];
rz(-2.5448749) q[2];
sx q[2];
rz(-2.2133568) q[2];
rz(-2.8440031) q[3];
sx q[3];
rz(-0.15407763) q[3];
sx q[3];
rz(0.84460622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069227844) q[0];
sx q[0];
rz(-2.897825) q[0];
sx q[0];
rz(0.099076554) q[0];
rz(-2.15436) q[1];
sx q[1];
rz(-1.3039373) q[1];
sx q[1];
rz(2.6027021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19463704) q[0];
sx q[0];
rz(-1.4879003) q[0];
sx q[0];
rz(3.1350053) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.107223) q[2];
sx q[2];
rz(-3.0863783) q[2];
sx q[2];
rz(-2.5718073) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3618631) q[1];
sx q[1];
rz(-1.2113844) q[1];
sx q[1];
rz(2.5796109) q[1];
rz(1.2480368) q[3];
sx q[3];
rz(-0.071167067) q[3];
sx q[3];
rz(-0.038480345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.146356) q[2];
sx q[2];
rz(-3.1260999) q[2];
sx q[2];
rz(0.30919477) q[2];
rz(0.63979465) q[3];
sx q[3];
rz(-0.00034172405) q[3];
sx q[3];
rz(0.063808002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30534202) q[0];
sx q[0];
rz(-2.5449365) q[0];
sx q[0];
rz(0.054656595) q[0];
rz(1.1326185) q[1];
sx q[1];
rz(-1.9839958) q[1];
sx q[1];
rz(1.2861015) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0564467) q[0];
sx q[0];
rz(-2.9582199) q[0];
sx q[0];
rz(-1.8118748) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.040820597) q[2];
sx q[2];
rz(-1.6128745) q[2];
sx q[2];
rz(3.1067318) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.4234679) q[1];
sx q[1];
rz(-1.7959975) q[1];
sx q[1];
rz(-1.4783212) q[1];
x q[2];
rz(-0.07225424) q[3];
sx q[3];
rz(-2.0569306) q[3];
sx q[3];
rz(0.89943258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9201811) q[2];
sx q[2];
rz(-2.5765918) q[2];
sx q[2];
rz(2.0378225) q[2];
rz(1.53299) q[3];
sx q[3];
rz(-3.0975603) q[3];
sx q[3];
rz(-2.485763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.00103818) q[0];
sx q[0];
rz(-2.9619205) q[0];
sx q[0];
rz(0.0044862577) q[0];
rz(1.5890315) q[1];
sx q[1];
rz(-1.6922502) q[1];
sx q[1];
rz(-3.0844614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2343952) q[0];
sx q[0];
rz(-2.937754) q[0];
sx q[0];
rz(-1.0782526) q[0];
rz(-1.7148027) q[2];
sx q[2];
rz(-1.3504354) q[2];
sx q[2];
rz(-1.2978467) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0100952) q[1];
sx q[1];
rz(-0.68298498) q[1];
sx q[1];
rz(-1.7630153) q[1];
x q[2];
rz(-0.83052333) q[3];
sx q[3];
rz(-1.4581513) q[3];
sx q[3];
rz(0.5121246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.39157465) q[2];
sx q[2];
rz(-3.1142758) q[2];
sx q[2];
rz(0.86891437) q[2];
rz(-2.1598375) q[3];
sx q[3];
rz(-0.029564094) q[3];
sx q[3];
rz(-2.6386007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045573087) q[0];
sx q[0];
rz(-1.4843142) q[0];
sx q[0];
rz(1.6577161) q[0];
rz(2.6847196) q[1];
sx q[1];
rz(-0.15468205) q[1];
sx q[1];
rz(-0.044943132) q[1];
rz(-0.76404608) q[2];
sx q[2];
rz(-1.7052787) q[2];
sx q[2];
rz(1.6324001) q[2];
rz(-1.7140688) q[3];
sx q[3];
rz(-0.91268703) q[3];
sx q[3];
rz(0.042162195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
