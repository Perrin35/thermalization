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
rz(1.60496032238007) q[0];
sx q[0];
rz(3.30242382188375) q[0];
sx q[0];
rz(9.81605801581546) q[0];
rz(-0.043035089969635) q[1];
sx q[1];
rz(4.76247623761232) q[1];
sx q[1];
rz(12.2279920339505) q[1];
cx q[1],q[0];
rz(1.79263770580292) q[0];
sx q[0];
rz(2.75119388301904) q[0];
sx q[0];
rz(7.0428185224454) q[0];
rz(-3.37023687362671) q[2];
sx q[2];
rz(6.50324407418305) q[2];
sx q[2];
rz(9.11724049448177) q[2];
cx q[2],q[1];
rz(-0.307928115129471) q[1];
sx q[1];
rz(5.37650814850862) q[1];
sx q[1];
rz(14.3796696424405) q[1];
rz(1.15427994728088) q[3];
sx q[3];
rz(4.81182626088197) q[3];
sx q[3];
rz(12.3276934385221) q[3];
cx q[3],q[2];
rz(0.368386209011078) q[2];
sx q[2];
rz(2.16190043290193) q[2];
sx q[2];
rz(8.366776084892) q[2];
rz(3.03927850723267) q[3];
sx q[3];
rz(4.75917843182618) q[3];
sx q[3];
rz(10.8250875234525) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.82329869270325) q[0];
sx q[0];
rz(2.38663909037644) q[0];
sx q[0];
rz(9.63867228328391) q[0];
rz(-1.33380436897278) q[1];
sx q[1];
rz(5.78294077714021) q[1];
sx q[1];
rz(12.2776069402616) q[1];
cx q[1],q[0];
rz(3.78174376487732) q[0];
sx q[0];
rz(4.08759883244569) q[0];
sx q[0];
rz(10.7826133727948) q[0];
rz(2.2367217540741) q[2];
sx q[2];
rz(2.18734321196611) q[2];
sx q[2];
rz(8.21118078231021) q[2];
cx q[2],q[1];
rz(-0.604644775390625) q[1];
sx q[1];
rz(4.68274012406404) q[1];
sx q[1];
rz(9.64561528562709) q[1];
rz(0.337776899337769) q[3];
sx q[3];
rz(-1.15725168387359) q[3];
sx q[3];
rz(12.4687862157743) q[3];
cx q[3],q[2];
rz(2.52209615707397) q[2];
sx q[2];
rz(2.32138261397416) q[2];
sx q[2];
rz(8.14043221472904) q[2];
rz(1.8340824842453) q[3];
sx q[3];
rz(4.45922807057435) q[3];
sx q[3];
rz(8.7262936592023) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.4035325050354) q[0];
sx q[0];
rz(4.20811775525148) q[0];
sx q[0];
rz(8.57948312758609) q[0];
rz(5.37760591506958) q[1];
sx q[1];
rz(3.74653610785539) q[1];
sx q[1];
rz(14.5303134679715) q[1];
cx q[1],q[0];
rz(2.98634099960327) q[0];
sx q[0];
rz(6.38029208977754) q[0];
sx q[0];
rz(8.06067214011356) q[0];
rz(4.18880796432495) q[2];
sx q[2];
rz(5.19959464867646) q[2];
sx q[2];
rz(9.95244613885089) q[2];
cx q[2],q[1];
rz(2.66034531593323) q[1];
sx q[1];
rz(1.77789095242555) q[1];
sx q[1];
rz(5.29325530528232) q[1];
rz(-0.645117819309235) q[3];
sx q[3];
rz(4.84019461472566) q[3];
sx q[3];
rz(10.0367125034253) q[3];
cx q[3],q[2];
rz(-2.52510690689087) q[2];
sx q[2];
rz(5.33898130257661) q[2];
sx q[2];
rz(9.29416823982402) q[2];
rz(1.7836000919342) q[3];
sx q[3];
rz(5.2744678576761) q[3];
sx q[3];
rz(11.2783679723661) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.65508770942688) q[0];
sx q[0];
rz(6.01885786850984) q[0];
sx q[0];
rz(12.9804720640103) q[0];
rz(0.865302383899689) q[1];
sx q[1];
rz(5.04620674450929) q[1];
sx q[1];
rz(10.0280374646108) q[1];
cx q[1],q[0];
rz(0.719724595546722) q[0];
sx q[0];
rz(2.65289211471612) q[0];
sx q[0];
rz(10.0575451016347) q[0];
rz(-3.07049536705017) q[2];
sx q[2];
rz(4.36120870907838) q[2];
sx q[2];
rz(13.0629119634549) q[2];
cx q[2],q[1];
rz(-1.32842600345612) q[1];
sx q[1];
rz(6.63405028183992) q[1];
sx q[1];
rz(11.1207837819974) q[1];
rz(2.78496813774109) q[3];
sx q[3];
rz(6.61544433434541) q[3];
sx q[3];
rz(10.7548303365628) q[3];
cx q[3],q[2];
rz(-3.66571021080017) q[2];
sx q[2];
rz(4.67393544514711) q[2];
sx q[2];
rz(11.9069626092832) q[2];
rz(3.00903296470642) q[3];
sx q[3];
rz(3.68826016982133) q[3];
sx q[3];
rz(11.8621434926908) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.346831262111664) q[0];
sx q[0];
rz(2.18573877413804) q[0];
sx q[0];
rz(9.15551007389232) q[0];
rz(1.5003833770752) q[1];
sx q[1];
rz(5.00413742859895) q[1];
sx q[1];
rz(11.4868590593259) q[1];
cx q[1],q[0];
rz(-3.52004337310791) q[0];
sx q[0];
rz(5.35200897057588) q[0];
sx q[0];
rz(12.7903716325681) q[0];
rz(-5.49039316177368) q[2];
sx q[2];
rz(2.07688203652436) q[2];
sx q[2];
rz(12.7702173948209) q[2];
cx q[2],q[1];
rz(-2.77918457984924) q[1];
sx q[1];
rz(5.19147911866242) q[1];
sx q[1];
rz(7.4790693283002) q[1];
rz(1.88873505592346) q[3];
sx q[3];
rz(1.77221813996369) q[3];
sx q[3];
rz(7.96655497550174) q[3];
cx q[3],q[2];
rz(-2.42148041725159) q[2];
sx q[2];
rz(3.72188261349732) q[2];
sx q[2];
rz(8.55454329251453) q[2];
rz(-1.80577111244202) q[3];
sx q[3];
rz(4.47709790070588) q[3];
sx q[3];
rz(8.74905732869312) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.133833527565) q[0];
sx q[0];
rz(3.1133911238336) q[0];
sx q[0];
rz(8.81355378626987) q[0];
rz(2.408034324646) q[1];
sx q[1];
rz(4.61240688164765) q[1];
sx q[1];
rz(9.80808893441364) q[1];
cx q[1],q[0];
rz(4.22833490371704) q[0];
sx q[0];
rz(5.40463844140107) q[0];
sx q[0];
rz(12.9306239843289) q[0];
rz(1.36938881874084) q[2];
sx q[2];
rz(3.97918620904023) q[2];
sx q[2];
rz(12.2016265153806) q[2];
cx q[2],q[1];
rz(5.34138059616089) q[1];
sx q[1];
rz(-0.933362332982473) q[1];
sx q[1];
rz(14.2965073347013) q[1];
rz(-0.844489395618439) q[3];
sx q[3];
rz(-1.25251135031646) q[3];
sx q[3];
rz(10.1431494712751) q[3];
cx q[3],q[2];
rz(-0.402062594890594) q[2];
sx q[2];
rz(2.58670875628526) q[2];
sx q[2];
rz(14.6950425863187) q[2];
rz(0.597645103931427) q[3];
sx q[3];
rz(4.87419977982576) q[3];
sx q[3];
rz(10.3507387995641) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.05858540534973) q[0];
sx q[0];
rz(3.57077220280702) q[0];
sx q[0];
rz(12.6023032426755) q[0];
rz(1.2220470905304) q[1];
sx q[1];
rz(3.81759593089158) q[1];
sx q[1];
rz(9.67284978031322) q[1];
cx q[1],q[0];
rz(0.00294686760753393) q[0];
sx q[0];
rz(-0.322975007695607) q[0];
sx q[0];
rz(10.0720227718274) q[0];
rz(1.08011341094971) q[2];
sx q[2];
rz(4.58731153805787) q[2];
sx q[2];
rz(9.4408072180967) q[2];
cx q[2],q[1];
rz(-2.38630199432373) q[1];
sx q[1];
rz(-0.404980031651906) q[1];
sx q[1];
rz(9.85013861059352) q[1];
rz(-0.00255360291339457) q[3];
sx q[3];
rz(4.81732693512971) q[3];
sx q[3];
rz(7.88464925288364) q[3];
cx q[3],q[2];
rz(1.0849357843399) q[2];
sx q[2];
rz(4.20129481156404) q[2];
sx q[2];
rz(5.75108573435947) q[2];
rz(-2.68712401390076) q[3];
sx q[3];
rz(1.54585090477998) q[3];
sx q[3];
rz(7.5772150516431) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.716595113277435) q[0];
sx q[0];
rz(2.09185615380342) q[0];
sx q[0];
rz(9.6649671703498) q[0];
rz(-0.605130851268768) q[1];
sx q[1];
rz(4.93541458447511) q[1];
sx q[1];
rz(10.0710999131124) q[1];
cx q[1],q[0];
rz(-1.54515647888184) q[0];
sx q[0];
rz(5.0120725949579) q[0];
sx q[0];
rz(11.8233737707059) q[0];
rz(0.0971664041280746) q[2];
sx q[2];
rz(2.77879134018952) q[2];
sx q[2];
rz(10.8150391340177) q[2];
cx q[2],q[1];
rz(-2.53902149200439) q[1];
sx q[1];
rz(2.63843098481233) q[1];
sx q[1];
rz(14.5723528623502) q[1];
rz(-1.78704917430878) q[3];
sx q[3];
rz(4.06825486023957) q[3];
sx q[3];
rz(8.60426745413944) q[3];
cx q[3],q[2];
rz(1.02950823307037) q[2];
sx q[2];
rz(4.78954020340974) q[2];
sx q[2];
rz(6.20150372981235) q[2];
rz(-0.847808480262756) q[3];
sx q[3];
rz(3.55085730751092) q[3];
sx q[3];
rz(12.3434066533963) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.341437011957169) q[0];
sx q[0];
rz(5.43849483330781) q[0];
sx q[0];
rz(10.3164592146794) q[0];
rz(-1.91012299060822) q[1];
sx q[1];
rz(0.488580854731151) q[1];
sx q[1];
rz(11.9565171956937) q[1];
cx q[1],q[0];
rz(1.11407399177551) q[0];
sx q[0];
rz(2.43540123303468) q[0];
sx q[0];
rz(8.66943333148166) q[0];
rz(1.58534038066864) q[2];
sx q[2];
rz(2.35682800610597) q[2];
sx q[2];
rz(9.97047779559299) q[2];
cx q[2],q[1];
rz(5.61468029022217) q[1];
sx q[1];
rz(5.44372216065461) q[1];
sx q[1];
rz(9.82281304001018) q[1];
rz(3.09503126144409) q[3];
sx q[3];
rz(3.48704239924485) q[3];
sx q[3];
rz(9.96811196803256) q[3];
cx q[3],q[2];
rz(1.86511981487274) q[2];
sx q[2];
rz(4.10888591607148) q[2];
sx q[2];
rz(12.836990094177) q[2];
rz(0.545177698135376) q[3];
sx q[3];
rz(4.95535746415193) q[3];
sx q[3];
rz(9.67294981180831) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.644286274909973) q[0];
sx q[0];
rz(4.44932106335694) q[0];
sx q[0];
rz(9.23828843831226) q[0];
rz(1.20937740802765) q[1];
sx q[1];
rz(4.49702897866304) q[1];
sx q[1];
rz(9.4051269184719) q[1];
cx q[1],q[0];
rz(-1.60158967971802) q[0];
sx q[0];
rz(2.18129095633561) q[0];
sx q[0];
rz(8.3592396736066) q[0];
rz(-3.41054272651672) q[2];
sx q[2];
rz(7.6017166693979) q[2];
sx q[2];
rz(8.69084272383853) q[2];
cx q[2],q[1];
rz(-0.841485917568207) q[1];
sx q[1];
rz(1.72386136849458) q[1];
sx q[1];
rz(11.351226425163) q[1];
rz(0.912274241447449) q[3];
sx q[3];
rz(2.26114747126634) q[3];
sx q[3];
rz(9.46645853518649) q[3];
cx q[3],q[2];
rz(0.432021170854568) q[2];
sx q[2];
rz(4.09705731471116) q[2];
sx q[2];
rz(7.8355842590253) q[2];
rz(0.109136953949928) q[3];
sx q[3];
rz(1.86926320393617) q[3];
sx q[3];
rz(9.148419058315) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.84211373329163) q[0];
sx q[0];
rz(3.50956899126107) q[0];
sx q[0];
rz(9.16440317629977) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(4.65705156326294) q[1];
sx q[1];
rz(2.58255586226518) q[1];
sx q[1];
rz(11.1751856565396) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-1.86024069786072) q[2];
sx q[2];
rz(4.70406297047669) q[2];
sx q[2];
rz(9.85425711273357) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.58692514896393) q[3];
sx q[3];
rz(1.95165148575837) q[3];
sx q[3];
rz(14.2330312490384) q[3];
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
