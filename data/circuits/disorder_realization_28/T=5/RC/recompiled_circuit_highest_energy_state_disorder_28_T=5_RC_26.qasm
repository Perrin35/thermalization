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
rz(-2.2313843) q[0];
sx q[0];
rz(3.892133) q[0];
sx q[0];
rz(12.005796) q[0];
rz(1.9380467) q[1];
sx q[1];
rz(-2.7653341) q[1];
sx q[1];
rz(2.9484152) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5995418) q[0];
sx q[0];
rz(-0.99630574) q[0];
sx q[0];
rz(-1.0380787) q[0];
rz(-pi) q[1];
rz(-1.4273663) q[2];
sx q[2];
rz(-0.71038112) q[2];
sx q[2];
rz(1.3262891) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3456125) q[1];
sx q[1];
rz(-2.4761845) q[1];
sx q[1];
rz(-1.7980952) q[1];
rz(-pi) q[2];
rz(-2.6704392) q[3];
sx q[3];
rz(-0.36590016) q[3];
sx q[3];
rz(0.19972502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91118497) q[2];
sx q[2];
rz(-2.513803) q[2];
sx q[2];
rz(0.5578624) q[2];
rz(1.2281536) q[3];
sx q[3];
rz(-1.4598673) q[3];
sx q[3];
rz(3.0465928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.8519583) q[0];
sx q[0];
rz(-1.292922) q[0];
sx q[0];
rz(1.5648382) q[0];
rz(-1.9948888) q[1];
sx q[1];
rz(-1.7161938) q[1];
sx q[1];
rz(-1.0316521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1250336) q[0];
sx q[0];
rz(-1.3368589) q[0];
sx q[0];
rz(2.1159322) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41298683) q[2];
sx q[2];
rz(-1.3527762) q[2];
sx q[2];
rz(-0.08139164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90742753) q[1];
sx q[1];
rz(-1.9112559) q[1];
sx q[1];
rz(-0.84788604) q[1];
rz(-0.35801631) q[3];
sx q[3];
rz(-2.6081134) q[3];
sx q[3];
rz(1.8915576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32776323) q[2];
sx q[2];
rz(-0.29080614) q[2];
sx q[2];
rz(-1.4168868) q[2];
rz(0.82301569) q[3];
sx q[3];
rz(-1.5282642) q[3];
sx q[3];
rz(-0.64457646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7773892) q[0];
sx q[0];
rz(-1.1193898) q[0];
sx q[0];
rz(-0.35225824) q[0];
rz(2.7525355) q[1];
sx q[1];
rz(-1.16951) q[1];
sx q[1];
rz(-0.22638098) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8181747) q[0];
sx q[0];
rz(-1.3238088) q[0];
sx q[0];
rz(3.0564708) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1445072) q[2];
sx q[2];
rz(-1.8663532) q[2];
sx q[2];
rz(1.3546582) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76629868) q[1];
sx q[1];
rz(-2.0336478) q[1];
sx q[1];
rz(2.0867324) q[1];
rz(-pi) q[2];
rz(2.705615) q[3];
sx q[3];
rz(-1.6395901) q[3];
sx q[3];
rz(1.1747557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9809197) q[2];
sx q[2];
rz(-1.0272762) q[2];
sx q[2];
rz(2.7799907) q[2];
rz(-2.8412039) q[3];
sx q[3];
rz(-1.3114248) q[3];
sx q[3];
rz(-0.96295199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2196197) q[0];
sx q[0];
rz(-2.0491056) q[0];
sx q[0];
rz(-0.12119448) q[0];
rz(-1.9425862) q[1];
sx q[1];
rz(-1.0795178) q[1];
sx q[1];
rz(1.8669063) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5124487) q[0];
sx q[0];
rz(-1.3487909) q[0];
sx q[0];
rz(-1.52284) q[0];
rz(0.077353296) q[2];
sx q[2];
rz(-2.935545) q[2];
sx q[2];
rz(-1.2543343) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.53654799) q[1];
sx q[1];
rz(-2.138461) q[1];
sx q[1];
rz(-1.3090538) q[1];
x q[2];
rz(-1.9363112) q[3];
sx q[3];
rz(-1.4344707) q[3];
sx q[3];
rz(-2.8773235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.092209665) q[2];
sx q[2];
rz(-1.694724) q[2];
sx q[2];
rz(1.764074) q[2];
rz(-2.0130017) q[3];
sx q[3];
rz(-2.1994574) q[3];
sx q[3];
rz(-2.3142864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89151299) q[0];
sx q[0];
rz(-0.87123195) q[0];
sx q[0];
rz(-1.0846035) q[0];
rz(2.4705823) q[1];
sx q[1];
rz(-1.5120466) q[1];
sx q[1];
rz(3.0674518) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0947572) q[0];
sx q[0];
rz(-2.1291385) q[0];
sx q[0];
rz(-1.2461016) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2535278) q[2];
sx q[2];
rz(-1.2937577) q[2];
sx q[2];
rz(-2.679897) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2977777) q[1];
sx q[1];
rz(-0.98688302) q[1];
sx q[1];
rz(2.816114) q[1];
x q[2];
rz(-1.1187333) q[3];
sx q[3];
rz(-2.8376355) q[3];
sx q[3];
rz(0.077465103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1159749) q[2];
sx q[2];
rz(-0.74883777) q[2];
sx q[2];
rz(-0.96308127) q[2];
rz(-1.7668308) q[3];
sx q[3];
rz(-2.0120967) q[3];
sx q[3];
rz(-0.53205427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2087723) q[0];
sx q[0];
rz(-2.8498579) q[0];
sx q[0];
rz(-1.8836841) q[0];
rz(-0.84016291) q[1];
sx q[1];
rz(-1.1779307) q[1];
sx q[1];
rz(0.69200969) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1245354) q[0];
sx q[0];
rz(-1.3003674) q[0];
sx q[0];
rz(-1.9084683) q[0];
rz(2.3756723) q[2];
sx q[2];
rz(-1.4408852) q[2];
sx q[2];
rz(-0.15405642) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6868478) q[1];
sx q[1];
rz(-2.1917289) q[1];
sx q[1];
rz(-2.5140155) q[1];
rz(0.63661544) q[3];
sx q[3];
rz(-2.0354384) q[3];
sx q[3];
rz(2.1615504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.96364) q[2];
sx q[2];
rz(-2.0452363) q[2];
sx q[2];
rz(-0.5160416) q[2];
rz(1.3456723) q[3];
sx q[3];
rz(-0.44107744) q[3];
sx q[3];
rz(1.0015944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8733785) q[0];
sx q[0];
rz(-2.8554947) q[0];
sx q[0];
rz(-1.0299261) q[0];
rz(-2.4870807) q[1];
sx q[1];
rz(-1.8067358) q[1];
sx q[1];
rz(-0.17394224) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2692341) q[0];
sx q[0];
rz(-2.1790811) q[0];
sx q[0];
rz(2.48003) q[0];
rz(-pi) q[1];
rz(-0.51474606) q[2];
sx q[2];
rz(-2.2651849) q[2];
sx q[2];
rz(-1.0389164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0566114) q[1];
sx q[1];
rz(-0.95889303) q[1];
sx q[1];
rz(-0.33501825) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6967322) q[3];
sx q[3];
rz(-2.5813817) q[3];
sx q[3];
rz(-1.0154533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30764636) q[2];
sx q[2];
rz(-0.53457326) q[2];
sx q[2];
rz(-1.3391116) q[2];
rz(0.80859679) q[3];
sx q[3];
rz(-1.9924889) q[3];
sx q[3];
rz(-1.8554525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22063743) q[0];
sx q[0];
rz(-2.0957102) q[0];
sx q[0];
rz(1.165423) q[0];
rz(2.9479345) q[1];
sx q[1];
rz(-0.8332738) q[1];
sx q[1];
rz(-1.9906893) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26955596) q[0];
sx q[0];
rz(-1.317476) q[0];
sx q[0];
rz(0.66590683) q[0];
rz(-pi) q[1];
rz(-1.3940349) q[2];
sx q[2];
rz(-0.77406787) q[2];
sx q[2];
rz(1.7947527) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85291687) q[1];
sx q[1];
rz(-1.1467458) q[1];
sx q[1];
rz(2.1217594) q[1];
rz(-0.70827534) q[3];
sx q[3];
rz(-2.357956) q[3];
sx q[3];
rz(-2.2275138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.2151486) q[2];
sx q[2];
rz(-1.0215267) q[2];
sx q[2];
rz(1.8683757) q[2];
rz(2.6269954) q[3];
sx q[3];
rz(-0.319258) q[3];
sx q[3];
rz(0.74015051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99106818) q[0];
sx q[0];
rz(-3.0758698) q[0];
sx q[0];
rz(2.7142628) q[0];
rz(-1.7006251) q[1];
sx q[1];
rz(-1.7040375) q[1];
sx q[1];
rz(0.89868054) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60481614) q[0];
sx q[0];
rz(-1.1425848) q[0];
sx q[0];
rz(0.34128071) q[0];
rz(-pi) q[1];
rz(-0.91941339) q[2];
sx q[2];
rz(-0.69910895) q[2];
sx q[2];
rz(-1.0369911) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9153122) q[1];
sx q[1];
rz(-1.8050005) q[1];
sx q[1];
rz(2.9620902) q[1];
rz(2.4472404) q[3];
sx q[3];
rz(-2.4789841) q[3];
sx q[3];
rz(-0.41678762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.906189) q[2];
sx q[2];
rz(-2.0960505) q[2];
sx q[2];
rz(-1.40353) q[2];
rz(2.4710726) q[3];
sx q[3];
rz(-1.5454005) q[3];
sx q[3];
rz(1.3129354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-0.32187605) q[0];
sx q[0];
rz(-0.96915594) q[0];
sx q[0];
rz(0.72625351) q[0];
rz(-0.67604524) q[1];
sx q[1];
rz(-1.36422) q[1];
sx q[1];
rz(2.3647251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99771927) q[0];
sx q[0];
rz(-1.5251499) q[0];
sx q[0];
rz(-0.93046988) q[0];
x q[1];
rz(-1.5258342) q[2];
sx q[2];
rz(-1.7167712) q[2];
sx q[2];
rz(-0.0073429664) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.54884262) q[1];
sx q[1];
rz(-1.934086) q[1];
sx q[1];
rz(0.32925683) q[1];
x q[2];
rz(-0.59263521) q[3];
sx q[3];
rz(-1.0909972) q[3];
sx q[3];
rz(2.0247839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1399416) q[2];
sx q[2];
rz(-0.17639128) q[2];
sx q[2];
rz(1.6528992) q[2];
rz(2.0148924) q[3];
sx q[3];
rz(-1.9727581) q[3];
sx q[3];
rz(2.6039629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61459944) q[0];
sx q[0];
rz(-0.64170964) q[0];
sx q[0];
rz(1.5668305) q[0];
rz(-0.78631403) q[1];
sx q[1];
rz(-0.17335261) q[1];
sx q[1];
rz(-0.39549624) q[1];
rz(2.4695071) q[2];
sx q[2];
rz(-1.0955878) q[2];
sx q[2];
rz(0.53447117) q[2];
rz(-0.24675225) q[3];
sx q[3];
rz(-2.3507531) q[3];
sx q[3];
rz(1.7274461) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
