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
rz(4.9024138) q[1];
sx q[1];
rz(9.496357) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78808331) q[0];
sx q[0];
rz(-3.1147163) q[0];
sx q[0];
rz(0.36969097) q[0];
rz(-pi) q[1];
rz(-1.3355005) q[2];
sx q[2];
rz(-1.1898196) q[2];
sx q[2];
rz(-1.345696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6759807) q[1];
sx q[1];
rz(-3.1221924) q[1];
sx q[1];
rz(-1.8288724) q[1];
rz(1.8911701) q[3];
sx q[3];
rz(-1.8719058) q[3];
sx q[3];
rz(-1.4442577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2804395) q[2];
sx q[2];
rz(-2.4344378) q[2];
sx q[2];
rz(-2.6748778) q[2];
rz(-2.6672582) q[3];
sx q[3];
rz(-3.1204087) q[3];
sx q[3];
rz(3.1202313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83653432) q[0];
sx q[0];
rz(-0.49764043) q[0];
sx q[0];
rz(-3.1217788) q[0];
rz(1.548798) q[1];
sx q[1];
rz(-2.9217547) q[1];
sx q[1];
rz(1.4733431) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9112355) q[0];
sx q[0];
rz(-0.94382554) q[0];
sx q[0];
rz(0.66642739) q[0];
x q[1];
rz(-0.85313932) q[2];
sx q[2];
rz(-0.38333508) q[2];
sx q[2];
rz(2.9789973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8791318) q[1];
sx q[1];
rz(-1.7647129) q[1];
sx q[1];
rz(-0.078207774) q[1];
x q[2];
rz(0.83062828) q[3];
sx q[3];
rz(-0.98139709) q[3];
sx q[3];
rz(2.7089861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3026498) q[2];
sx q[2];
rz(-2.5431716) q[2];
sx q[2];
rz(-1.2897276) q[2];
rz(-1.9047811) q[3];
sx q[3];
rz(-2.808282) q[3];
sx q[3];
rz(0.70241565) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1692093) q[0];
sx q[0];
rz(-1.1595668) q[0];
sx q[0];
rz(1.5706536) q[0];
rz(1.4717357) q[1];
sx q[1];
rz(-1.5262628) q[1];
sx q[1];
rz(-0.45438802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4226625) q[0];
sx q[0];
rz(-0.61394982) q[0];
sx q[0];
rz(-1.3387247) q[0];
x q[1];
rz(1.387792) q[2];
sx q[2];
rz(-2.9972074) q[2];
sx q[2];
rz(2.3600277) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66167605) q[1];
sx q[1];
rz(-3.0212086) q[1];
sx q[1];
rz(1.3811587) q[1];
rz(-pi) q[2];
rz(-0.64735712) q[3];
sx q[3];
rz(-1.0203927) q[3];
sx q[3];
rz(3.0542706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4001974) q[2];
sx q[2];
rz(-0.016308451) q[2];
sx q[2];
rz(-0.34747094) q[2];
rz(-2.6853284) q[3];
sx q[3];
rz(-3.1268692) q[3];
sx q[3];
rz(0.99304503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149813) q[0];
sx q[0];
rz(-1.240629) q[0];
sx q[0];
rz(1.4062784) q[0];
rz(-0.44334385) q[1];
sx q[1];
rz(-1.025238) q[1];
sx q[1];
rz(-1.5700856) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3298108) q[0];
sx q[0];
rz(-2.499281) q[0];
sx q[0];
rz(-0.091600939) q[0];
rz(-pi) q[1];
rz(-2.3894044) q[2];
sx q[2];
rz(-3.0471932) q[2];
sx q[2];
rz(-2.029325) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0896596) q[1];
sx q[1];
rz(-1.5610236) q[1];
sx q[1];
rz(1.8490514) q[1];
x q[2];
rz(1.1394386) q[3];
sx q[3];
rz(-1.316081) q[3];
sx q[3];
rz(-1.6540119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37322474) q[2];
sx q[2];
rz(-2.7478605) q[2];
sx q[2];
rz(0.079252871) q[2];
rz(1.9521889) q[3];
sx q[3];
rz(-1.7704084) q[3];
sx q[3];
rz(-1.6325379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.70334148) q[0];
sx q[0];
rz(-0.51333135) q[0];
sx q[0];
rz(-2.3355423) q[0];
rz(-0.84000677) q[1];
sx q[1];
rz(-0.012931074) q[1];
sx q[1];
rz(-2.3654225) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989202) q[0];
sx q[0];
rz(-1.4399043) q[0];
sx q[0];
rz(-2.5517273) q[0];
rz(-pi) q[1];
rz(3.0018535) q[2];
sx q[2];
rz(-0.010541803) q[2];
sx q[2];
rz(0.14103061) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8766206) q[1];
sx q[1];
rz(-1.7062467) q[1];
sx q[1];
rz(-1.5828062) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25960323) q[3];
sx q[3];
rz(-1.4624551) q[3];
sx q[3];
rz(1.3110127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1157896) q[2];
sx q[2];
rz(-1.5805406) q[2];
sx q[2];
rz(-0.72581327) q[2];
rz(-2.9535661) q[3];
sx q[3];
rz(-0.060421061) q[3];
sx q[3];
rz(-2.4316725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96227729) q[0];
sx q[0];
rz(-0.55914068) q[0];
sx q[0];
rz(0.54139262) q[0];
rz(-2.9560282) q[1];
sx q[1];
rz(-1.5927529) q[1];
sx q[1];
rz(-3.013179) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1379515) q[0];
sx q[0];
rz(-1.9310474) q[0];
sx q[0];
rz(1.2482578) q[0];
x q[1];
rz(-1.4494684) q[2];
sx q[2];
rz(-1.575496) q[2];
sx q[2];
rz(3.1374861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6278319) q[1];
sx q[1];
rz(-2.7134088) q[1];
sx q[1];
rz(2.9437888) q[1];
rz(1.6261934) q[3];
sx q[3];
rz(-0.18871103) q[3];
sx q[3];
rz(-0.95079899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7554756) q[2];
sx q[2];
rz(-0.057689276) q[2];
sx q[2];
rz(2.3003787) q[2];
rz(2.9403213) q[3];
sx q[3];
rz(-1.5320675) q[3];
sx q[3];
rz(-2.8819528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5614618) q[0];
sx q[0];
rz(-2.353297) q[0];
sx q[0];
rz(-1.5607675) q[0];
rz(0.67281094) q[1];
sx q[1];
rz(-1.7377868) q[1];
sx q[1];
rz(-3.1127081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7772693) q[0];
sx q[0];
rz(-1.8176846) q[0];
sx q[0];
rz(0.29233513) q[0];
rz(0.47841623) q[2];
sx q[2];
rz(-1.3923651) q[2];
sx q[2];
rz(-2.894424) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.72994443) q[1];
sx q[1];
rz(-1.8081417) q[1];
sx q[1];
rz(-1.1470331) q[1];
rz(-2.6050909) q[3];
sx q[3];
rz(-2.287545) q[3];
sx q[3];
rz(1.052618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2334571) q[2];
sx q[2];
rz(-2.5448749) q[2];
sx q[2];
rz(2.2133568) q[2];
rz(0.2975896) q[3];
sx q[3];
rz(-2.987515) q[3];
sx q[3];
rz(-0.84460622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069227844) q[0];
sx q[0];
rz(-0.2437676) q[0];
sx q[0];
rz(-3.0425161) q[0];
rz(-2.15436) q[1];
sx q[1];
rz(-1.3039373) q[1];
sx q[1];
rz(2.6027021) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7659788) q[0];
sx q[0];
rz(-1.5642316) q[0];
sx q[0];
rz(1.6536941) q[0];
rz(2.0343696) q[2];
sx q[2];
rz(-3.0863783) q[2];
sx q[2];
rz(0.56978536) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7797295) q[1];
sx q[1];
rz(-1.9302082) q[1];
sx q[1];
rz(-0.56198175) q[1];
rz(-pi) q[2];
rz(1.8935558) q[3];
sx q[3];
rz(-3.0704256) q[3];
sx q[3];
rz(3.1031123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.146356) q[2];
sx q[2];
rz(-0.015492798) q[2];
sx q[2];
rz(2.8323979) q[2];
rz(-2.501798) q[3];
sx q[3];
rz(-0.00034172405) q[3];
sx q[3];
rz(-3.0777847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30534202) q[0];
sx q[0];
rz(-0.59665614) q[0];
sx q[0];
rz(-0.054656595) q[0];
rz(-2.0089741) q[1];
sx q[1];
rz(-1.1575969) q[1];
sx q[1];
rz(-1.2861015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7228381) q[0];
sx q[0];
rz(-1.5272473) q[0];
sx q[0];
rz(1.7489793) q[0];
x q[1];
rz(-3.1007721) q[2];
sx q[2];
rz(-1.6128745) q[2];
sx q[2];
rz(0.034860858) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1117592) q[1];
sx q[1];
rz(-0.24315029) q[1];
sx q[1];
rz(-2.7584226) q[1];
x q[2];
rz(1.70657) q[3];
sx q[3];
rz(-2.6505436) q[3];
sx q[3];
rz(0.74573475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9201811) q[2];
sx q[2];
rz(-2.5765918) q[2];
sx q[2];
rz(1.1037702) q[2];
rz(1.53299) q[3];
sx q[3];
rz(-3.0975603) q[3];
sx q[3];
rz(-2.485763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.00103818) q[0];
sx q[0];
rz(-2.9619205) q[0];
sx q[0];
rz(0.0044862577) q[0];
rz(1.5525612) q[1];
sx q[1];
rz(-1.6922502) q[1];
sx q[1];
rz(3.0844614) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.321314) q[0];
sx q[0];
rz(-1.6666659) q[0];
sx q[0];
rz(1.3906327) q[0];
x q[1];
rz(0.56978858) q[2];
sx q[2];
rz(-0.26260153) q[2];
sx q[2];
rz(-0.71209967) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.58913402) q[1];
sx q[1];
rz(-1.4499364) q[1];
sx q[1];
rz(-2.2446988) q[1];
rz(-pi) q[2];
rz(1.404625) q[3];
sx q[3];
rz(-0.7471841) q[3];
sx q[3];
rz(0.93624828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39157465) q[2];
sx q[2];
rz(-0.027316814) q[2];
sx q[2];
rz(-0.86891437) q[2];
rz(-0.9817552) q[3];
sx q[3];
rz(-3.1120286) q[3];
sx q[3];
rz(0.50299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.045573087) q[0];
sx q[0];
rz(-1.6572784) q[0];
sx q[0];
rz(-1.4838765) q[0];
rz(2.6847196) q[1];
sx q[1];
rz(-0.15468205) q[1];
sx q[1];
rz(-0.044943132) q[1];
rz(2.3775466) q[2];
sx q[2];
rz(-1.7052787) q[2];
sx q[2];
rz(1.6324001) q[2];
rz(0.66309912) q[3];
sx q[3];
rz(-1.4575921) q[3];
sx q[3];
rz(1.7009638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
