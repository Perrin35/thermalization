OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91035834) q[0];
sx q[0];
rz(-2.2725821) q[0];
sx q[0];
rz(2.056871) q[0];
rz(1.9864858) q[1];
sx q[1];
rz(-2.3218563) q[1];
sx q[1];
rz(0.91135946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706575) q[0];
sx q[0];
rz(-0.80659272) q[0];
sx q[0];
rz(0.33814221) q[0];
x q[1];
rz(0.16884825) q[2];
sx q[2];
rz(-2.422794) q[2];
sx q[2];
rz(1.079725) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.45381308) q[1];
sx q[1];
rz(-1.1247083) q[1];
sx q[1];
rz(-0.58961726) q[1];
rz(-1.7098268) q[3];
sx q[3];
rz(-2.7043501) q[3];
sx q[3];
rz(2.4831877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.800941) q[2];
sx q[2];
rz(-1.1110577) q[2];
sx q[2];
rz(-0.079455376) q[2];
rz(-0.60845145) q[3];
sx q[3];
rz(-1.918101) q[3];
sx q[3];
rz(-2.2505545) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4873753) q[0];
sx q[0];
rz(-2.2779164) q[0];
sx q[0];
rz(-1.3551711) q[0];
rz(-1.0379418) q[1];
sx q[1];
rz(-2.4619921) q[1];
sx q[1];
rz(-0.80744809) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6896967) q[0];
sx q[0];
rz(-2.1679255) q[0];
sx q[0];
rz(0.35291617) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1218929) q[2];
sx q[2];
rz(-2.4906213) q[2];
sx q[2];
rz(-0.62007444) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8958853) q[1];
sx q[1];
rz(-1.504532) q[1];
sx q[1];
rz(-2.9016728) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1065527) q[3];
sx q[3];
rz(-1.7037183) q[3];
sx q[3];
rz(-0.35655856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21343931) q[2];
sx q[2];
rz(-1.5779053) q[2];
sx q[2];
rz(-1.6820924) q[2];
rz(-0.45808074) q[3];
sx q[3];
rz(-0.57772485) q[3];
sx q[3];
rz(2.1817082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15993519) q[0];
sx q[0];
rz(-1.0502879) q[0];
sx q[0];
rz(-2.4666069) q[0];
rz(0.58810294) q[1];
sx q[1];
rz(-2.2876078) q[1];
sx q[1];
rz(-1.3311707) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25615293) q[0];
sx q[0];
rz(-2.0567499) q[0];
sx q[0];
rz(1.11973) q[0];
rz(-2.0841062) q[2];
sx q[2];
rz(-0.9710487) q[2];
sx q[2];
rz(0.63569234) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5471897) q[1];
sx q[1];
rz(-1.8246027) q[1];
sx q[1];
rz(2.3188792) q[1];
x q[2];
rz(1.0717594) q[3];
sx q[3];
rz(-1.3924358) q[3];
sx q[3];
rz(0.87770977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62446928) q[2];
sx q[2];
rz(-0.34999592) q[2];
sx q[2];
rz(2.9208753) q[2];
rz(-0.80495009) q[3];
sx q[3];
rz(-1.5996108) q[3];
sx q[3];
rz(-1.3360924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6014366) q[0];
sx q[0];
rz(-0.24208459) q[0];
sx q[0];
rz(0.618774) q[0];
rz(3.0586808) q[1];
sx q[1];
rz(-2.7989048) q[1];
sx q[1];
rz(1.6212911) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12099685) q[0];
sx q[0];
rz(-1.6669287) q[0];
sx q[0];
rz(-2.0970048) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57915202) q[2];
sx q[2];
rz(-1.7699827) q[2];
sx q[2];
rz(2.3483089) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1136696) q[1];
sx q[1];
rz(-1.4884559) q[1];
sx q[1];
rz(-2.0827977) q[1];
rz(-2.7079775) q[3];
sx q[3];
rz(-2.7325897) q[3];
sx q[3];
rz(-0.65164372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.80359047) q[2];
sx q[2];
rz(-3.0323196) q[2];
sx q[2];
rz(-2.9962311) q[2];
rz(-2.8694782) q[3];
sx q[3];
rz(-1.7348671) q[3];
sx q[3];
rz(-1.9947778) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014932545) q[0];
sx q[0];
rz(-1.8155875) q[0];
sx q[0];
rz(-3.0354011) q[0];
rz(0.95942489) q[1];
sx q[1];
rz(-2.5966849) q[1];
sx q[1];
rz(-2.9471961) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1379762) q[0];
sx q[0];
rz(-2.4038393) q[0];
sx q[0];
rz(-0.12807782) q[0];
rz(-pi) q[1];
rz(1.4559559) q[2];
sx q[2];
rz(-1.5659589) q[2];
sx q[2];
rz(0.23813914) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.186708) q[1];
sx q[1];
rz(-0.59046035) q[1];
sx q[1];
rz(-0.073530274) q[1];
rz(-pi) q[2];
rz(3.0644981) q[3];
sx q[3];
rz(-1.1655679) q[3];
sx q[3];
rz(0.73540686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.306119) q[2];
sx q[2];
rz(-1.4543616) q[2];
sx q[2];
rz(2.5731738) q[2];
rz(-3.0889555) q[3];
sx q[3];
rz(-1.6391552) q[3];
sx q[3];
rz(3.0047825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1133872) q[0];
sx q[0];
rz(-0.55332342) q[0];
sx q[0];
rz(2.9521039) q[0];
rz(-1.0264617) q[1];
sx q[1];
rz(-2.2920513) q[1];
sx q[1];
rz(-0.75526563) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.301195) q[0];
sx q[0];
rz(-0.9009255) q[0];
sx q[0];
rz(-2.2706881) q[0];
x q[1];
rz(2.3113475) q[2];
sx q[2];
rz(-0.70607042) q[2];
sx q[2];
rz(-0.65078562) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6615852) q[1];
sx q[1];
rz(-0.23218854) q[1];
sx q[1];
rz(2.3228374) q[1];
rz(-0.35474687) q[3];
sx q[3];
rz(-1.2257396) q[3];
sx q[3];
rz(-1.3416895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3834164) q[2];
sx q[2];
rz(-1.5864317) q[2];
sx q[2];
rz(1.4413393) q[2];
rz(1.3759184) q[3];
sx q[3];
rz(-2.223189) q[3];
sx q[3];
rz(-2.4413696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43948424) q[0];
sx q[0];
rz(-2.0937884) q[0];
sx q[0];
rz(-1.1442319) q[0];
rz(-2.8817835) q[1];
sx q[1];
rz(-0.83994284) q[1];
sx q[1];
rz(-0.98794404) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081962498) q[0];
sx q[0];
rz(-3.0407627) q[0];
sx q[0];
rz(-2.9743845) q[0];
x q[1];
rz(2.4056882) q[2];
sx q[2];
rz(-1.2469075) q[2];
sx q[2];
rz(1.0867556) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0212342) q[1];
sx q[1];
rz(-1.8250006) q[1];
sx q[1];
rz(-0.66616504) q[1];
x q[2];
rz(-2.1454303) q[3];
sx q[3];
rz(-1.8271433) q[3];
sx q[3];
rz(-0.72130313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7320431) q[2];
sx q[2];
rz(-1.5012375) q[2];
sx q[2];
rz(0.6130971) q[2];
rz(2.4260855) q[3];
sx q[3];
rz(-0.77176538) q[3];
sx q[3];
rz(0.36090052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9382984) q[0];
sx q[0];
rz(-0.84380117) q[0];
sx q[0];
rz(0.26813689) q[0];
rz(-2.9109491) q[1];
sx q[1];
rz(-1.5985039) q[1];
sx q[1];
rz(-3.1275829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21523968) q[0];
sx q[0];
rz(-0.82414675) q[0];
sx q[0];
rz(1.8570379) q[0];
rz(-1.2410937) q[2];
sx q[2];
rz(-0.46247855) q[2];
sx q[2];
rz(-1.3792737) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4625293) q[1];
sx q[1];
rz(-2.6025297) q[1];
sx q[1];
rz(-0.087058914) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1118891) q[3];
sx q[3];
rz(-1.4280768) q[3];
sx q[3];
rz(-2.5102455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9739428) q[2];
sx q[2];
rz(-1.8029982) q[2];
sx q[2];
rz(-2.8323284) q[2];
rz(-2.9850128) q[3];
sx q[3];
rz(-1.4045818) q[3];
sx q[3];
rz(-1.0369161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.88368791) q[0];
sx q[0];
rz(-2.4990999) q[0];
sx q[0];
rz(2.8908492) q[0];
rz(0.58654395) q[1];
sx q[1];
rz(-1.5778912) q[1];
sx q[1];
rz(-2.3768545) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4737807) q[0];
sx q[0];
rz(-0.43536738) q[0];
sx q[0];
rz(2.3171168) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3735869) q[2];
sx q[2];
rz(-0.80942488) q[2];
sx q[2];
rz(-2.1999733) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8426399) q[1];
sx q[1];
rz(-0.12315006) q[1];
sx q[1];
rz(0.075249057) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13651092) q[3];
sx q[3];
rz(-2.7595466) q[3];
sx q[3];
rz(2.7979224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6278729) q[2];
sx q[2];
rz(-0.79078117) q[2];
sx q[2];
rz(-2.9676843) q[2];
rz(-2.3993313) q[3];
sx q[3];
rz(-0.75209135) q[3];
sx q[3];
rz(-3.098587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-1.0088418) q[0];
sx q[0];
rz(-0.77021563) q[0];
sx q[0];
rz(-2.8736864) q[0];
rz(1.335089) q[1];
sx q[1];
rz(-2.2309512) q[1];
sx q[1];
rz(1.7693899) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89248449) q[0];
sx q[0];
rz(-1.3297446) q[0];
sx q[0];
rz(-1.1066828) q[0];
rz(2.82136) q[2];
sx q[2];
rz(-0.31117421) q[2];
sx q[2];
rz(2.8428889) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7497014) q[1];
sx q[1];
rz(-1.0969583) q[1];
sx q[1];
rz(2.5577765) q[1];
rz(-pi) q[2];
rz(-1.2413512) q[3];
sx q[3];
rz(-1.2678896) q[3];
sx q[3];
rz(-2.0249413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.043896) q[2];
sx q[2];
rz(-1.0363657) q[2];
sx q[2];
rz(1.3807266) q[2];
rz(2.0128287) q[3];
sx q[3];
rz(-0.94996101) q[3];
sx q[3];
rz(-1.5560163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632521) q[0];
sx q[0];
rz(-1.5338407) q[0];
sx q[0];
rz(-1.1080678) q[0];
rz(-0.68645984) q[1];
sx q[1];
rz(-1.3931128) q[1];
sx q[1];
rz(-1.211094) q[1];
rz(1.2800693) q[2];
sx q[2];
rz(-1.0940934) q[2];
sx q[2];
rz(1.4296503) q[2];
rz(-2.7298418) q[3];
sx q[3];
rz(-2.0709548) q[3];
sx q[3];
rz(-1.0695468) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
