OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.546917200088501) q[0];
sx q[0];
rz(4.14247992833192) q[0];
sx q[0];
rz(9.63717978297874) q[0];
rz(0.714959323406219) q[1];
sx q[1];
rz(3.92908045847947) q[1];
sx q[1];
rz(10.7062834262769) q[1];
cx q[1],q[0];
rz(0.403853833675385) q[0];
sx q[0];
rz(4.91325107415254) q[0];
sx q[0];
rz(9.07745454310581) q[0];
rz(-0.878269374370575) q[2];
sx q[2];
rz(1.5480969270044) q[2];
sx q[2];
rz(11.3412396669309) q[2];
cx q[2],q[1];
rz(1.23385298252106) q[1];
sx q[1];
rz(5.82360878785188) q[1];
sx q[1];
rz(12.9519409894864) q[1];
rz(0.12870566546917) q[3];
sx q[3];
rz(4.77500382264192) q[3];
sx q[3];
rz(8.96535820364162) q[3];
cx q[3],q[2];
rz(-0.32221257686615) q[2];
sx q[2];
rz(3.25654314656789) q[2];
sx q[2];
rz(10.9415130376737) q[2];
rz(-0.247622817754745) q[3];
sx q[3];
rz(4.61534169514711) q[3];
sx q[3];
rz(11.9108178377072) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.44419288635254) q[0];
sx q[0];
rz(4.80510917504365) q[0];
sx q[0];
rz(9.5810718446891) q[0];
rz(0.639317572116852) q[1];
sx q[1];
rz(5.27057877381379) q[1];
sx q[1];
rz(7.7720243692319) q[1];
cx q[1],q[0];
rz(-0.811936795711517) q[0];
sx q[0];
rz(2.56630793412263) q[0];
sx q[0];
rz(10.6999773740689) q[0];
rz(0.747295081615448) q[2];
sx q[2];
rz(4.03739330370957) q[2];
sx q[2];
rz(8.37908241748019) q[2];
cx q[2],q[1];
rz(1.38581514358521) q[1];
sx q[1];
rz(7.28725353081758) q[1];
sx q[1];
rz(10.3328842282216) q[1];
rz(1.35272467136383) q[3];
sx q[3];
rz(3.28146651585633) q[3];
sx q[3];
rz(10.5281810522) q[3];
cx q[3],q[2];
rz(-0.850410103797913) q[2];
sx q[2];
rz(3.24047107447917) q[2];
sx q[2];
rz(9.12806636690303) q[2];
rz(0.332427799701691) q[3];
sx q[3];
rz(4.7647022326761) q[3];
sx q[3];
rz(10.6586131811063) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.63337028026581) q[0];
sx q[0];
rz(4.40035870869691) q[0];
sx q[0];
rz(9.95509318112537) q[0];
rz(-0.5439532995224) q[1];
sx q[1];
rz(4.26302424271638) q[1];
sx q[1];
rz(10.6138099193494) q[1];
cx q[1],q[0];
rz(1.31845676898956) q[0];
sx q[0];
rz(1.41347769101197) q[0];
sx q[0];
rz(11.2952641010205) q[0];
rz(1.72789549827576) q[2];
sx q[2];
rz(3.39274040062959) q[2];
sx q[2];
rz(9.05412462948962) q[2];
cx q[2],q[1];
rz(-0.584547340869904) q[1];
sx q[1];
rz(4.18555227120454) q[1];
sx q[1];
rz(13.3396606206815) q[1];
rz(0.872394442558289) q[3];
sx q[3];
rz(1.93653705914552) q[3];
sx q[3];
rz(10.9114479780118) q[3];
cx q[3],q[2];
rz(2.23607802391052) q[2];
sx q[2];
rz(4.81017783482606) q[2];
sx q[2];
rz(10.676338171951) q[2];
rz(0.454231470823288) q[3];
sx q[3];
rz(2.70827037294442) q[3];
sx q[3];
rz(9.77991283535167) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.778047561645508) q[0];
sx q[0];
rz(3.45128193696076) q[0];
sx q[0];
rz(10.765719151489) q[0];
rz(0.606079339981079) q[1];
sx q[1];
rz(5.34416619141633) q[1];
sx q[1];
rz(10.6094282627027) q[1];
cx q[1],q[0];
rz(-0.236833244562149) q[0];
sx q[0];
rz(1.39688113530213) q[0];
sx q[0];
rz(8.07161507605716) q[0];
rz(0.72261905670166) q[2];
sx q[2];
rz(3.78531798918779) q[2];
sx q[2];
rz(11.0704236984174) q[2];
cx q[2],q[1];
rz(0.468900978565216) q[1];
sx q[1];
rz(5.68117085297639) q[1];
sx q[1];
rz(12.7734088659207) q[1];
rz(-0.366551458835602) q[3];
sx q[3];
rz(4.31281379063661) q[3];
sx q[3];
rz(10.6279361009519) q[3];
cx q[3],q[2];
rz(-0.233453154563904) q[2];
sx q[2];
rz(5.39737835724885) q[2];
sx q[2];
rz(10.0643416404645) q[2];
rz(2.28416395187378) q[3];
sx q[3];
rz(1.99684170086915) q[3];
sx q[3];
rz(8.93498367666408) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.630675911903381) q[0];
sx q[0];
rz(2.59987214406068) q[0];
sx q[0];
rz(12.0459770917813) q[0];
rz(-0.102945417165756) q[1];
sx q[1];
rz(4.19204035599763) q[1];
sx q[1];
rz(10.1535892248075) q[1];
cx q[1],q[0];
rz(1.67293524742126) q[0];
sx q[0];
rz(2.91420826514299) q[0];
sx q[0];
rz(10.546609735481) q[0];
rz(-0.183766335248947) q[2];
sx q[2];
rz(2.4090320785814) q[2];
sx q[2];
rz(11.0114384651105) q[2];
cx q[2],q[1];
rz(1.55069577693939) q[1];
sx q[1];
rz(4.50787785847718) q[1];
sx q[1];
rz(9.75414264797374) q[1];
rz(0.613522469997406) q[3];
sx q[3];
rz(4.25606349309022) q[3];
sx q[3];
rz(9.858605092756) q[3];
cx q[3],q[2];
rz(1.09854280948639) q[2];
sx q[2];
rz(3.50151208241517) q[2];
sx q[2];
rz(6.99289772509738) q[2];
rz(0.633061587810516) q[3];
sx q[3];
rz(4.28975466092164) q[3];
sx q[3];
rz(9.29306500255271) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.10280406475067) q[0];
sx q[0];
rz(6.18081274827058) q[0];
sx q[0];
rz(10.3084950804631) q[0];
rz(0.920641303062439) q[1];
sx q[1];
rz(5.22110787232453) q[1];
sx q[1];
rz(6.48083446025058) q[1];
cx q[1],q[0];
rz(0.550068557262421) q[0];
sx q[0];
rz(5.32704153855378) q[0];
sx q[0];
rz(10.2399065852086) q[0];
rz(1.46428239345551) q[2];
sx q[2];
rz(4.79134860833222) q[2];
sx q[2];
rz(11.8874191999356) q[2];
cx q[2],q[1];
rz(0.270539879798889) q[1];
sx q[1];
rz(3.59856102068956) q[1];
sx q[1];
rz(10.8988642454068) q[1];
rz(-0.686244606971741) q[3];
sx q[3];
rz(5.34587732155854) q[3];
sx q[3];
rz(10.7462657451551) q[3];
cx q[3],q[2];
rz(0.0278230421245098) q[2];
sx q[2];
rz(4.93967834313447) q[2];
sx q[2];
rz(12.6419749021451) q[2];
rz(1.48937010765076) q[3];
sx q[3];
rz(3.54106974800164) q[3];
sx q[3];
rz(5.45410511492893) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.419121146202087) q[0];
sx q[0];
rz(4.02706721623475) q[0];
sx q[0];
rz(9.4668510839264) q[0];
rz(1.17885899543762) q[1];
sx q[1];
rz(1.98316398461396) q[1];
sx q[1];
rz(11.396938776962) q[1];
cx q[1],q[0];
rz(1.62993156909943) q[0];
sx q[0];
rz(3.31941887934739) q[0];
sx q[0];
rz(9.81060779689952) q[0];
rz(-1.96096539497375) q[2];
sx q[2];
rz(4.97621944745118) q[2];
sx q[2];
rz(8.52216295003101) q[2];
cx q[2],q[1];
rz(0.840502560138702) q[1];
sx q[1];
rz(1.23553481896455) q[1];
sx q[1];
rz(10.4236793875615) q[1];
rz(-3.09939527511597) q[3];
sx q[3];
rz(1.95846668084199) q[3];
sx q[3];
rz(11.2736195087354) q[3];
cx q[3],q[2];
rz(-1.82326698303223) q[2];
sx q[2];
rz(2.02647581894929) q[2];
sx q[2];
rz(9.43926131948038) q[2];
rz(1.62183403968811) q[3];
sx q[3];
rz(4.43838420708711) q[3];
sx q[3];
rz(8.87159851788684) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.139772146940231) q[0];
sx q[0];
rz(2.79123899539048) q[0];
sx q[0];
rz(9.98712018727466) q[0];
rz(-1.43158042430878) q[1];
sx q[1];
rz(4.88619390328462) q[1];
sx q[1];
rz(9.88251284360095) q[1];
cx q[1],q[0];
rz(-0.916775286197662) q[0];
sx q[0];
rz(3.62528190215165) q[0];
sx q[0];
rz(9.24423857628509) q[0];
rz(0.696763932704926) q[2];
sx q[2];
rz(4.27120676835115) q[2];
sx q[2];
rz(11.6886863469998) q[2];
cx q[2],q[1];
rz(2.88439583778381) q[1];
sx q[1];
rz(5.44909468491609) q[1];
sx q[1];
rz(7.8346510887067) q[1];
rz(-0.551499128341675) q[3];
sx q[3];
rz(5.68443790276582) q[3];
sx q[3];
rz(12.1951133966367) q[3];
cx q[3],q[2];
rz(0.835950136184692) q[2];
sx q[2];
rz(2.20174107153947) q[2];
sx q[2];
rz(9.70731524228259) q[2];
rz(2.18411684036255) q[3];
sx q[3];
rz(4.53383198578889) q[3];
sx q[3];
rz(8.94603282808467) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.958650588989258) q[0];
sx q[0];
rz(3.82689002354676) q[0];
sx q[0];
rz(11.1640772581021) q[0];
rz(-0.699334025382996) q[1];
sx q[1];
rz(4.97990122635896) q[1];
sx q[1];
rz(11.3045629024427) q[1];
cx q[1],q[0];
rz(0.487550735473633) q[0];
sx q[0];
rz(5.76834359963471) q[0];
sx q[0];
rz(11.5661813974301) q[0];
rz(0.917557120323181) q[2];
sx q[2];
rz(4.878482015925) q[2];
sx q[2];
rz(8.03078815936252) q[2];
cx q[2],q[1];
rz(-3.03144335746765) q[1];
sx q[1];
rz(4.45984628994996) q[1];
sx q[1];
rz(14.8419885396878) q[1];
rz(-1.0412175655365) q[3];
sx q[3];
rz(4.53452447255189) q[3];
sx q[3];
rz(11.4247047662656) q[3];
cx q[3],q[2];
rz(0.703014731407166) q[2];
sx q[2];
rz(4.64546755154664) q[2];
sx q[2];
rz(10.8826267480771) q[2];
rz(0.666493773460388) q[3];
sx q[3];
rz(4.56946709950502) q[3];
sx q[3];
rz(10.1764958858411) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.895019710063934) q[0];
sx q[0];
rz(3.95197031100328) q[0];
sx q[0];
rz(10.5847247600476) q[0];
rz(-0.0541404224932194) q[1];
sx q[1];
rz(1.47779467900331) q[1];
sx q[1];
rz(10.4952180147092) q[1];
cx q[1],q[0];
rz(2.11500716209412) q[0];
sx q[0];
rz(2.38376489480073) q[0];
sx q[0];
rz(9.66859615444347) q[0];
rz(-0.0398633480072021) q[2];
sx q[2];
rz(3.9067800958925) q[2];
sx q[2];
rz(12.0553471803586) q[2];
cx q[2],q[1];
rz(1.14203190803528) q[1];
sx q[1];
rz(3.97241363127763) q[1];
sx q[1];
rz(5.9442405462186) q[1];
rz(0.959941625595093) q[3];
sx q[3];
rz(5.98354783852632) q[3];
sx q[3];
rz(9.933470940582) q[3];
cx q[3],q[2];
rz(3.91035437583923) q[2];
sx q[2];
rz(5.73412719567353) q[2];
sx q[2];
rz(10.9516194820325) q[2];
rz(2.56201434135437) q[3];
sx q[3];
rz(4.93635073502595) q[3];
sx q[3];
rz(10.0768841862599) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.02544605731964) q[0];
sx q[0];
rz(2.28852918942506) q[0];
sx q[0];
rz(9.37609501033231) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.00289940834045) q[1];
sx q[1];
rz(2.00787857373292) q[1];
sx q[1];
rz(7.3017627954404) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.0760956257581711) q[2];
sx q[2];
rz(4.66970887978608) q[2];
sx q[2];
rz(10.1572621226232) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.215307742357254) q[3];
sx q[3];
rz(2.15874412854249) q[3];
sx q[3];
rz(13.1004848241727) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
