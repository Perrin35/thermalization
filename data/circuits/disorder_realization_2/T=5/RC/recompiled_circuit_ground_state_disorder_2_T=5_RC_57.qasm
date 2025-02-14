OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.53192294) q[0];
sx q[0];
rz(1.6165531) q[0];
sx q[0];
rz(11.308148) q[0];
rz(7.8426709) q[1];
sx q[1];
rz(5.7962228) q[1];
sx q[1];
rz(16.143057) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5601016) q[0];
sx q[0];
rz(-1.035443) q[0];
sx q[0];
rz(-2.125432) q[0];
rz(-2.6039106) q[2];
sx q[2];
rz(-1.2721561) q[2];
sx q[2];
rz(-0.60288211) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.62578979) q[1];
sx q[1];
rz(-2.8608659) q[1];
sx q[1];
rz(1.8516385) q[1];
rz(-0.55032309) q[3];
sx q[3];
rz(-1.1512412) q[3];
sx q[3];
rz(1.5418095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9293999) q[2];
sx q[2];
rz(-1.8636899) q[2];
sx q[2];
rz(1.0431935) q[2];
rz(2.7405401) q[3];
sx q[3];
rz(-1.4570823) q[3];
sx q[3];
rz(2.5636087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.9143518) q[0];
sx q[0];
rz(-0.14509097) q[0];
sx q[0];
rz(2.5703854) q[0];
rz(0.97194833) q[1];
sx q[1];
rz(-0.74058878) q[1];
sx q[1];
rz(-2.0170225) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6475424) q[0];
sx q[0];
rz(-2.4801773) q[0];
sx q[0];
rz(1.1986092) q[0];
x q[1];
rz(-2.6827181) q[2];
sx q[2];
rz(-0.70415184) q[2];
sx q[2];
rz(-1.6353232) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31370762) q[1];
sx q[1];
rz(-1.4769625) q[1];
sx q[1];
rz(2.3052144) q[1];
rz(-pi) q[2];
x q[2];
rz(1.589682) q[3];
sx q[3];
rz(-1.369993) q[3];
sx q[3];
rz(0.41057472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9198415) q[2];
sx q[2];
rz(-1.5566748) q[2];
sx q[2];
rz(0.43608967) q[2];
rz(0.74887216) q[3];
sx q[3];
rz(-0.80513969) q[3];
sx q[3];
rz(-2.3124636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1083168) q[0];
sx q[0];
rz(-1.746614) q[0];
sx q[0];
rz(3.0535611) q[0];
rz(-0.85405675) q[1];
sx q[1];
rz(-0.66103926) q[1];
sx q[1];
rz(1.0850272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4861761) q[0];
sx q[0];
rz(-1.2987505) q[0];
sx q[0];
rz(0.06975091) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7872115) q[2];
sx q[2];
rz(-1.4206352) q[2];
sx q[2];
rz(0.79126287) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4702556) q[1];
sx q[1];
rz(-0.47154676) q[1];
sx q[1];
rz(-1.9370473) q[1];
rz(2.964193) q[3];
sx q[3];
rz(-2.7466603) q[3];
sx q[3];
rz(-1.653914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4714841) q[2];
sx q[2];
rz(-1.4205616) q[2];
sx q[2];
rz(-0.74618435) q[2];
rz(0.21607312) q[3];
sx q[3];
rz(-1.8862628) q[3];
sx q[3];
rz(-2.9358673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7699319) q[0];
sx q[0];
rz(-0.07337229) q[0];
sx q[0];
rz(-1.7727456) q[0];
rz(2.5313077) q[1];
sx q[1];
rz(-1.7758324) q[1];
sx q[1];
rz(2.761421) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73715009) q[0];
sx q[0];
rz(-1.2594873) q[0];
sx q[0];
rz(0.4510757) q[0];
rz(2.6908003) q[2];
sx q[2];
rz(-1.3935699) q[2];
sx q[2];
rz(1.4361567) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7138765) q[1];
sx q[1];
rz(-0.52449924) q[1];
sx q[1];
rz(1.8668411) q[1];
rz(-pi) q[2];
x q[2];
rz(0.065297619) q[3];
sx q[3];
rz(-1.1698128) q[3];
sx q[3];
rz(2.9915265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2569106) q[2];
sx q[2];
rz(-2.1941049) q[2];
sx q[2];
rz(2.103629) q[2];
rz(-2.7773618) q[3];
sx q[3];
rz(-0.7887775) q[3];
sx q[3];
rz(-2.4558333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6895741) q[0];
sx q[0];
rz(-1.7114534) q[0];
sx q[0];
rz(2.3110287) q[0];
rz(3.0914302) q[1];
sx q[1];
rz(-3.0144033) q[1];
sx q[1];
rz(0.62279472) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1019331) q[0];
sx q[0];
rz(-2.9103511) q[0];
sx q[0];
rz(0.34666611) q[0];
rz(-pi) q[1];
rz(0.41983126) q[2];
sx q[2];
rz(-1.3908885) q[2];
sx q[2];
rz(-0.77984389) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0633358) q[1];
sx q[1];
rz(-1.2163269) q[1];
sx q[1];
rz(-0.88651231) q[1];
rz(-pi) q[2];
rz(0.022449724) q[3];
sx q[3];
rz(-1.5707309) q[3];
sx q[3];
rz(1.8409468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.043776) q[2];
sx q[2];
rz(-2.706683) q[2];
sx q[2];
rz(-0.80769509) q[2];
rz(-1.5378599) q[3];
sx q[3];
rz(-1.6350428) q[3];
sx q[3];
rz(0.22411552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1767126) q[0];
sx q[0];
rz(-1.3310615) q[0];
sx q[0];
rz(-1.4759195) q[0];
rz(-1.2337947) q[1];
sx q[1];
rz(-1.4101615) q[1];
sx q[1];
rz(2.9880611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3447405) q[0];
sx q[0];
rz(-1.5670876) q[0];
sx q[0];
rz(0.10964236) q[0];
x q[1];
rz(-0.45120181) q[2];
sx q[2];
rz(-2.0236733) q[2];
sx q[2];
rz(1.42746) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8922032) q[1];
sx q[1];
rz(-1.0703328) q[1];
sx q[1];
rz(1.0618147) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12795398) q[3];
sx q[3];
rz(-0.64264402) q[3];
sx q[3];
rz(1.2817739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4393356) q[2];
sx q[2];
rz(-2.3372237) q[2];
sx q[2];
rz(0.55629936) q[2];
rz(-3.124681) q[3];
sx q[3];
rz(-2.1664679) q[3];
sx q[3];
rz(-2.4771966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9282114) q[0];
sx q[0];
rz(-0.072755486) q[0];
sx q[0];
rz(2.4640006) q[0];
rz(1.3202336) q[1];
sx q[1];
rz(-1.0878539) q[1];
sx q[1];
rz(0.56603146) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91221936) q[0];
sx q[0];
rz(-2.006833) q[0];
sx q[0];
rz(-1.1993809) q[0];
rz(-3.1268372) q[2];
sx q[2];
rz(-1.5185673) q[2];
sx q[2];
rz(1.9697389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85105079) q[1];
sx q[1];
rz(-0.7194964) q[1];
sx q[1];
rz(-1.3403203) q[1];
rz(-pi) q[2];
rz(1.612408) q[3];
sx q[3];
rz(-1.0806314) q[3];
sx q[3];
rz(0.13652882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90938202) q[2];
sx q[2];
rz(-0.38572329) q[2];
sx q[2];
rz(2.5717226) q[2];
rz(2.0991523) q[3];
sx q[3];
rz(-1.9614377) q[3];
sx q[3];
rz(2.1377783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37782272) q[0];
sx q[0];
rz(-0.42709392) q[0];
sx q[0];
rz(2.1183993) q[0];
rz(-0.91776735) q[1];
sx q[1];
rz(-0.75420165) q[1];
sx q[1];
rz(0.86407152) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9167332) q[0];
sx q[0];
rz(-1.3274487) q[0];
sx q[0];
rz(0.12638522) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1064498) q[2];
sx q[2];
rz(-1.9111782) q[2];
sx q[2];
rz(-0.92591531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9814216) q[1];
sx q[1];
rz(-2.2653901) q[1];
sx q[1];
rz(0.92215529) q[1];
rz(-pi) q[2];
rz(-1.9404961) q[3];
sx q[3];
rz(-2.2636578) q[3];
sx q[3];
rz(0.38398352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0495905) q[2];
sx q[2];
rz(-1.2668173) q[2];
sx q[2];
rz(0.76466307) q[2];
rz(-2.2406254) q[3];
sx q[3];
rz(-2.4668906) q[3];
sx q[3];
rz(1.0990134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6242999) q[0];
sx q[0];
rz(-2.7462672) q[0];
sx q[0];
rz(-2.1584391) q[0];
rz(-0.82955018) q[1];
sx q[1];
rz(-1.364565) q[1];
sx q[1];
rz(0.57873908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5924517) q[0];
sx q[0];
rz(-1.5835174) q[0];
sx q[0];
rz(1.1770269) q[0];
rz(2.1095244) q[2];
sx q[2];
rz(-2.0499637) q[2];
sx q[2];
rz(0.62512809) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14565258) q[1];
sx q[1];
rz(-2.4697587) q[1];
sx q[1];
rz(0.89222676) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7962977) q[3];
sx q[3];
rz(-1.6052947) q[3];
sx q[3];
rz(-1.2277227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4911554) q[2];
sx q[2];
rz(-1.8530242) q[2];
sx q[2];
rz(-0.11162652) q[2];
rz(-2.1857183) q[3];
sx q[3];
rz(-2.1501232) q[3];
sx q[3];
rz(1.8808232) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2535506) q[0];
sx q[0];
rz(-2.084806) q[0];
sx q[0];
rz(3.1274617) q[0];
rz(-2.0290532) q[1];
sx q[1];
rz(-2.2281149) q[1];
sx q[1];
rz(1.6003476) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6131564) q[0];
sx q[0];
rz(-2.9455958) q[0];
sx q[0];
rz(-2.4264748) q[0];
x q[1];
rz(-0.12596031) q[2];
sx q[2];
rz(-0.93011127) q[2];
sx q[2];
rz(1.0116742) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4804876) q[1];
sx q[1];
rz(-2.7723162) q[1];
sx q[1];
rz(2.9758456) q[1];
rz(-pi) q[2];
rz(-2.9716206) q[3];
sx q[3];
rz(-1.408443) q[3];
sx q[3];
rz(1.8708097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9666226) q[2];
sx q[2];
rz(-1.9113767) q[2];
sx q[2];
rz(0.54565412) q[2];
rz(3.0555861) q[3];
sx q[3];
rz(-1.6274446) q[3];
sx q[3];
rz(2.8470794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1957112) q[0];
sx q[0];
rz(-2.0989037) q[0];
sx q[0];
rz(1.209191) q[0];
rz(-2.53881) q[1];
sx q[1];
rz(-2.5310015) q[1];
sx q[1];
rz(-2.3878154) q[1];
rz(0.88861619) q[2];
sx q[2];
rz(-1.4063514) q[2];
sx q[2];
rz(0.92264639) q[2];
rz(-0.62282563) q[3];
sx q[3];
rz(-1.6591743) q[3];
sx q[3];
rz(1.7521973) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
