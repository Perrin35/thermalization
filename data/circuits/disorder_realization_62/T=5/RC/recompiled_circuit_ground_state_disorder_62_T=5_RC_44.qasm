OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5247076) q[0];
sx q[0];
rz(-2.4029713) q[0];
sx q[0];
rz(-0.11697669) q[0];
rz(0.95247954) q[1];
sx q[1];
rz(-1.6713961) q[1];
sx q[1];
rz(-2.267946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3461733) q[0];
sx q[0];
rz(-1.367698) q[0];
sx q[0];
rz(0.83195306) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0181276) q[2];
sx q[2];
rz(-1.2510994) q[2];
sx q[2];
rz(-1.5716219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4651771) q[1];
sx q[1];
rz(-0.55182528) q[1];
sx q[1];
rz(1.4613749) q[1];
x q[2];
rz(1.9756894) q[3];
sx q[3];
rz(-1.8244074) q[3];
sx q[3];
rz(-0.028279956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.17254193) q[2];
sx q[2];
rz(-1.3323063) q[2];
sx q[2];
rz(1.8633899) q[2];
rz(1.5365907) q[3];
sx q[3];
rz(-1.9032109) q[3];
sx q[3];
rz(0.31755519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69739598) q[0];
sx q[0];
rz(-0.2187271) q[0];
sx q[0];
rz(0.87795192) q[0];
rz(2.3776993) q[1];
sx q[1];
rz(-2.0593675) q[1];
sx q[1];
rz(1.1306995) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3618523) q[0];
sx q[0];
rz(-1.6310143) q[0];
sx q[0];
rz(-1.4152933) q[0];
rz(-1.892852) q[2];
sx q[2];
rz(-2.6832504) q[2];
sx q[2];
rz(0.94142585) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2511828) q[1];
sx q[1];
rz(-1.5580388) q[1];
sx q[1];
rz(-2.3240672) q[1];
rz(1.6783707) q[3];
sx q[3];
rz(-0.94753303) q[3];
sx q[3];
rz(-3.0855765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6905602) q[2];
sx q[2];
rz(-1.6310383) q[2];
sx q[2];
rz(-0.65563273) q[2];
rz(-1.9111274) q[3];
sx q[3];
rz(-0.96223193) q[3];
sx q[3];
rz(-0.83278304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(1.7448298) q[0];
sx q[0];
rz(-1.1711045) q[0];
sx q[0];
rz(-0.50286621) q[0];
rz(-0.14547959) q[1];
sx q[1];
rz(-1.611404) q[1];
sx q[1];
rz(-1.9810289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.202464) q[0];
sx q[0];
rz(-1.2697497) q[0];
sx q[0];
rz(2.8168764) q[0];
rz(-2.7941452) q[2];
sx q[2];
rz(-0.37003741) q[2];
sx q[2];
rz(-1.7269788) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7532901) q[1];
sx q[1];
rz(-0.75725746) q[1];
sx q[1];
rz(2.1427335) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3480556) q[3];
sx q[3];
rz(-1.8249434) q[3];
sx q[3];
rz(2.4055907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9460556) q[2];
sx q[2];
rz(-1.035752) q[2];
sx q[2];
rz(-2.3717086) q[2];
rz(2.5464673) q[3];
sx q[3];
rz(-0.49134058) q[3];
sx q[3];
rz(-2.6814521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50764099) q[0];
sx q[0];
rz(-1.0469629) q[0];
sx q[0];
rz(-1.7858343) q[0];
rz(2.5178364) q[1];
sx q[1];
rz(-1.6056986) q[1];
sx q[1];
rz(-2.6502868) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2441933) q[0];
sx q[0];
rz(-1.4174875) q[0];
sx q[0];
rz(0.53863948) q[0];
x q[1];
rz(2.5042438) q[2];
sx q[2];
rz(-0.34301153) q[2];
sx q[2];
rz(-2.1528139) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92190336) q[1];
sx q[1];
rz(-1.2619201) q[1];
sx q[1];
rz(-0.83654735) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66670615) q[3];
sx q[3];
rz(-1.9940064) q[3];
sx q[3];
rz(-2.4413787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9451311) q[2];
sx q[2];
rz(-1.5245695) q[2];
sx q[2];
rz(2.5932236) q[2];
rz(-1.8114926) q[3];
sx q[3];
rz(-2.8434704) q[3];
sx q[3];
rz(-0.33838457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2746628) q[0];
sx q[0];
rz(-2.6882956) q[0];
sx q[0];
rz(-1.0484265) q[0];
rz(3.0323845) q[1];
sx q[1];
rz(-1.5914773) q[1];
sx q[1];
rz(1.106326) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8620548) q[0];
sx q[0];
rz(-1.3509271) q[0];
sx q[0];
rz(3.0310275) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8064617) q[2];
sx q[2];
rz(-2.2673528) q[2];
sx q[2];
rz(2.1256688) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8815461) q[1];
sx q[1];
rz(-1.4235919) q[1];
sx q[1];
rz(-1.24415) q[1];
rz(-pi) q[2];
rz(-2.2436687) q[3];
sx q[3];
rz(-2.7677288) q[3];
sx q[3];
rz(-1.3488692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76935524) q[2];
sx q[2];
rz(-1.2662338) q[2];
sx q[2];
rz(-2.9388169) q[2];
rz(-0.77962223) q[3];
sx q[3];
rz(-1.0896261) q[3];
sx q[3];
rz(-2.7705418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.62696537) q[0];
sx q[0];
rz(-0.045128673) q[0];
sx q[0];
rz(-2.2623862) q[0];
rz(-0.59576398) q[1];
sx q[1];
rz(-0.87409449) q[1];
sx q[1];
rz(1.2557868) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43508732) q[0];
sx q[0];
rz(-1.6235949) q[0];
sx q[0];
rz(-3.0493124) q[0];
x q[1];
rz(-2.6379616) q[2];
sx q[2];
rz(-0.63490311) q[2];
sx q[2];
rz(2.0201928) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3869898) q[1];
sx q[1];
rz(-2.5856254) q[1];
sx q[1];
rz(2.703859) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34481315) q[3];
sx q[3];
rz(-2.4969144) q[3];
sx q[3];
rz(0.0078474069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6229728) q[2];
sx q[2];
rz(-2.8039248) q[2];
sx q[2];
rz(-1.6597623) q[2];
rz(1.4431813) q[3];
sx q[3];
rz(-1.7549763) q[3];
sx q[3];
rz(-1.1186918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5959394) q[0];
sx q[0];
rz(-2.312199) q[0];
sx q[0];
rz(-2.948092) q[0];
rz(-0.045348383) q[1];
sx q[1];
rz(-1.9769316) q[1];
sx q[1];
rz(2.4605816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7398241) q[0];
sx q[0];
rz(-1.6575749) q[0];
sx q[0];
rz(-2.9784188) q[0];
rz(-2.0991285) q[2];
sx q[2];
rz(-1.7914504) q[2];
sx q[2];
rz(2.200732) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0630156) q[1];
sx q[1];
rz(-1.583644) q[1];
sx q[1];
rz(-0.0058369667) q[1];
rz(-0.38661031) q[3];
sx q[3];
rz(-1.6865239) q[3];
sx q[3];
rz(-0.9936617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8650774) q[2];
sx q[2];
rz(-2.4962208) q[2];
sx q[2];
rz(1.3481677) q[2];
rz(1.3639785) q[3];
sx q[3];
rz(-1.8230702) q[3];
sx q[3];
rz(-1.1905504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6337223) q[0];
sx q[0];
rz(-0.81955925) q[0];
sx q[0];
rz(-0.67935294) q[0];
rz(0.23718111) q[1];
sx q[1];
rz(-0.91865426) q[1];
sx q[1];
rz(2.0749345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54618597) q[0];
sx q[0];
rz(-2.0810938) q[0];
sx q[0];
rz(1.575281) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96644281) q[2];
sx q[2];
rz(-1.4928515) q[2];
sx q[2];
rz(1.9029531) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.28344) q[1];
sx q[1];
rz(-1.6951218) q[1];
sx q[1];
rz(1.0470301) q[1];
x q[2];
rz(-1.1447843) q[3];
sx q[3];
rz(-1.8717812) q[3];
sx q[3];
rz(2.3097484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25881585) q[2];
sx q[2];
rz(-1.9100185) q[2];
sx q[2];
rz(2.7302177) q[2];
rz(-1.8607633) q[3];
sx q[3];
rz(-2.4954093) q[3];
sx q[3];
rz(2.083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6537052) q[0];
sx q[0];
rz(-1.7937086) q[0];
sx q[0];
rz(-2.2109798) q[0];
rz(2.6764684) q[1];
sx q[1];
rz(-2.5587176) q[1];
sx q[1];
rz(-1.7238269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8450369) q[0];
sx q[0];
rz(-0.64528685) q[0];
sx q[0];
rz(0.40420239) q[0];
rz(-pi) q[1];
rz(0.82724656) q[2];
sx q[2];
rz(-0.8349903) q[2];
sx q[2];
rz(-0.5448676) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5807996) q[1];
sx q[1];
rz(-1.8617814) q[1];
sx q[1];
rz(-0.35620226) q[1];
x q[2];
rz(3.0259589) q[3];
sx q[3];
rz(-0.99471417) q[3];
sx q[3];
rz(-2.2779358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.96665367) q[2];
sx q[2];
rz(-0.83028364) q[2];
sx q[2];
rz(-1.2154328) q[2];
rz(-0.71839607) q[3];
sx q[3];
rz(-1.6738439) q[3];
sx q[3];
rz(-2.4992656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65828085) q[0];
sx q[0];
rz(-0.42228666) q[0];
sx q[0];
rz(-1.8102113) q[0];
rz(0.70679682) q[1];
sx q[1];
rz(-1.4247318) q[1];
sx q[1];
rz(-1.4250863) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1757705) q[0];
sx q[0];
rz(-0.81908161) q[0];
sx q[0];
rz(0.62721647) q[0];
rz(-pi) q[1];
rz(-0.064878929) q[2];
sx q[2];
rz(-0.91514313) q[2];
sx q[2];
rz(-0.82198373) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6725965) q[1];
sx q[1];
rz(-1.1332268) q[1];
sx q[1];
rz(-0.5794748) q[1];
rz(2.0277068) q[3];
sx q[3];
rz(-0.82172365) q[3];
sx q[3];
rz(0.81629074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66296545) q[2];
sx q[2];
rz(-0.68296432) q[2];
sx q[2];
rz(-1.9662201) q[2];
rz(-3.1183682) q[3];
sx q[3];
rz(-0.5718137) q[3];
sx q[3];
rz(0.2636675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57446734) q[0];
sx q[0];
rz(-2.1463285) q[0];
sx q[0];
rz(1.2405201) q[0];
rz(-0.65470882) q[1];
sx q[1];
rz(-1.2798825) q[1];
sx q[1];
rz(1.9725694) q[1];
rz(-0.11266795) q[2];
sx q[2];
rz(-1.4341528) q[2];
sx q[2];
rz(-1.2918593) q[2];
rz(-0.8169218) q[3];
sx q[3];
rz(-2.2754442) q[3];
sx q[3];
rz(1.2572921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
