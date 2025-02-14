OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8213788) q[0];
sx q[0];
rz(-2.3795405) q[0];
sx q[0];
rz(2.2574545) q[0];
rz(2.60131) q[1];
sx q[1];
rz(-1.5457431) q[1];
sx q[1];
rz(0.63313142) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0057982) q[0];
sx q[0];
rz(-1.457419) q[0];
sx q[0];
rz(-2.9762639) q[0];
rz(-pi) q[1];
rz(2.6548926) q[2];
sx q[2];
rz(-1.4632157) q[2];
sx q[2];
rz(2.8963367) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.84304436) q[1];
sx q[1];
rz(-1.6982276) q[1];
sx q[1];
rz(-0.28020466) q[1];
rz(-pi) q[2];
rz(-2.5930674) q[3];
sx q[3];
rz(-1.9270634) q[3];
sx q[3];
rz(-1.6690205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4151609) q[2];
sx q[2];
rz(-1.8981322) q[2];
sx q[2];
rz(2.4746056) q[2];
rz(0.64217448) q[3];
sx q[3];
rz(-2.8447076) q[3];
sx q[3];
rz(1.4750397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17382962) q[0];
sx q[0];
rz(-0.9742631) q[0];
sx q[0];
rz(0.94671384) q[0];
rz(0.058454839) q[1];
sx q[1];
rz(-2.7775601) q[1];
sx q[1];
rz(0.91190016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0094441) q[0];
sx q[0];
rz(-1.5329038) q[0];
sx q[0];
rz(-1.1542341) q[0];
rz(-pi) q[1];
x q[1];
rz(0.070621892) q[2];
sx q[2];
rz(-2.2433519) q[2];
sx q[2];
rz(0.20820752) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.81296125) q[1];
sx q[1];
rz(-1.8133303) q[1];
sx q[1];
rz(-0.46142908) q[1];
rz(-pi) q[2];
rz(-0.1095406) q[3];
sx q[3];
rz(-0.44407223) q[3];
sx q[3];
rz(-0.20694297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6902206) q[2];
sx q[2];
rz(-2.8896832) q[2];
sx q[2];
rz(0.33736649) q[2];
rz(0.99006027) q[3];
sx q[3];
rz(-1.3288574) q[3];
sx q[3];
rz(1.2216074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5628691) q[0];
sx q[0];
rz(-2.8674485) q[0];
sx q[0];
rz(1.3124056) q[0];
rz(-2.9263272) q[1];
sx q[1];
rz(-1.5039517) q[1];
sx q[1];
rz(0.85830918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90412437) q[0];
sx q[0];
rz(-1.6209813) q[0];
sx q[0];
rz(-1.3048569) q[0];
x q[1];
rz(-2.1444001) q[2];
sx q[2];
rz(-2.7421087) q[2];
sx q[2];
rz(-1.4956719) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2006113) q[1];
sx q[1];
rz(-1.1340144) q[1];
sx q[1];
rz(1.3882367) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86238548) q[3];
sx q[3];
rz(-0.38194712) q[3];
sx q[3];
rz(-0.89513429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12576088) q[2];
sx q[2];
rz(-0.55700004) q[2];
sx q[2];
rz(-1.0090984) q[2];
rz(0.5674181) q[3];
sx q[3];
rz(-1.6332473) q[3];
sx q[3];
rz(1.5198038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0030982) q[0];
sx q[0];
rz(-2.4479285) q[0];
sx q[0];
rz(0.19288119) q[0];
rz(-2.0774972) q[1];
sx q[1];
rz(-2.8654983) q[1];
sx q[1];
rz(2.1853133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67502695) q[0];
sx q[0];
rz(-2.1888632) q[0];
sx q[0];
rz(-1.7873554) q[0];
x q[1];
rz(-0.85727875) q[2];
sx q[2];
rz(-0.801221) q[2];
sx q[2];
rz(2.4080133) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5259486) q[1];
sx q[1];
rz(-0.61405776) q[1];
sx q[1];
rz(-2.581852) q[1];
rz(-pi) q[2];
rz(-0.90518454) q[3];
sx q[3];
rz(-1.815261) q[3];
sx q[3];
rz(0.28206952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8320147) q[2];
sx q[2];
rz(-0.39065209) q[2];
sx q[2];
rz(-2.2008956) q[2];
rz(-0.37788033) q[3];
sx q[3];
rz(-0.89972073) q[3];
sx q[3];
rz(-2.5009724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2864895) q[0];
sx q[0];
rz(-1.6096977) q[0];
sx q[0];
rz(1.8163649) q[0];
rz(2.8638966) q[1];
sx q[1];
rz(-2.1546021) q[1];
sx q[1];
rz(1.206254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2640058) q[0];
sx q[0];
rz(-0.20795272) q[0];
sx q[0];
rz(-1.1780187) q[0];
rz(-pi) q[1];
rz(1.9812859) q[2];
sx q[2];
rz(-2.1417981) q[2];
sx q[2];
rz(-2.0620907) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6305651) q[1];
sx q[1];
rz(-1.2113809) q[1];
sx q[1];
rz(0.014822249) q[1];
rz(-2.8998781) q[3];
sx q[3];
rz(-1.7346584) q[3];
sx q[3];
rz(-3.0047447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.18579927) q[2];
sx q[2];
rz(-2.9904371) q[2];
sx q[2];
rz(1.7929057) q[2];
rz(-1.0079314) q[3];
sx q[3];
rz(-1.4069087) q[3];
sx q[3];
rz(-3.1048043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0808829) q[0];
sx q[0];
rz(-3.1302629) q[0];
sx q[0];
rz(2.9445373) q[0];
rz(-1.9673678) q[1];
sx q[1];
rz(-0.86001992) q[1];
sx q[1];
rz(2.7161782) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2046609) q[0];
sx q[0];
rz(-2.6789894) q[0];
sx q[0];
rz(3.0499056) q[0];
x q[1];
rz(-1.4324911) q[2];
sx q[2];
rz(-0.30685654) q[2];
sx q[2];
rz(2.8371161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.08845932) q[1];
sx q[1];
rz(-1.3091374) q[1];
sx q[1];
rz(0.43105521) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86023902) q[3];
sx q[3];
rz(-1.0061534) q[3];
sx q[3];
rz(-2.3535239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60874385) q[2];
sx q[2];
rz(-0.11561919) q[2];
sx q[2];
rz(-0.70333424) q[2];
rz(0.15130791) q[3];
sx q[3];
rz(-1.1727099) q[3];
sx q[3];
rz(0.54733384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5587826) q[0];
sx q[0];
rz(-0.8081688) q[0];
sx q[0];
rz(2.1852344) q[0];
rz(-2.3725678) q[1];
sx q[1];
rz(-1.290657) q[1];
sx q[1];
rz(-2.1020611) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35339402) q[0];
sx q[0];
rz(-2.0457342) q[0];
sx q[0];
rz(0.026508412) q[0];
rz(-pi) q[1];
rz(0.65724397) q[2];
sx q[2];
rz(-2.1604156) q[2];
sx q[2];
rz(0.60962617) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7733147) q[1];
sx q[1];
rz(-1.8551143) q[1];
sx q[1];
rz(1.8892688) q[1];
rz(1.9045181) q[3];
sx q[3];
rz(-1.7705373) q[3];
sx q[3];
rz(-2.7488696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2018373) q[2];
sx q[2];
rz(-0.73362437) q[2];
sx q[2];
rz(-1.6377829) q[2];
rz(-0.87320295) q[3];
sx q[3];
rz(-2.2538897) q[3];
sx q[3];
rz(2.8100815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37209508) q[0];
sx q[0];
rz(-2.7517419) q[0];
sx q[0];
rz(1.1497585) q[0];
rz(0.86504966) q[1];
sx q[1];
rz(-1.9416092) q[1];
sx q[1];
rz(-1.7080151) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5652611) q[0];
sx q[0];
rz(-2.1339685) q[0];
sx q[0];
rz(-0.43656008) q[0];
x q[1];
rz(-1.9351129) q[2];
sx q[2];
rz(-1.120036) q[2];
sx q[2];
rz(0.70012602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3506602) q[1];
sx q[1];
rz(-2.6106842) q[1];
sx q[1];
rz(-1.3965308) q[1];
rz(-2.0074437) q[3];
sx q[3];
rz(-0.45387196) q[3];
sx q[3];
rz(-1.3144823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4992497) q[2];
sx q[2];
rz(-1.6982102) q[2];
sx q[2];
rz(0.7439822) q[2];
rz(2.3104525) q[3];
sx q[3];
rz(-1.7404514) q[3];
sx q[3];
rz(0.35541117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1019854) q[0];
sx q[0];
rz(-2.2672741) q[0];
sx q[0];
rz(2.4264477) q[0];
rz(-1.2932418) q[1];
sx q[1];
rz(-1.5512385) q[1];
sx q[1];
rz(0.23983461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3837357) q[0];
sx q[0];
rz(-0.78690398) q[0];
sx q[0];
rz(-0.73946865) q[0];
rz(-pi) q[1];
rz(1.5101899) q[2];
sx q[2];
rz(-2.1327634) q[2];
sx q[2];
rz(-2.2417703) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.71069392) q[1];
sx q[1];
rz(-1.7393117) q[1];
sx q[1];
rz(-2.2546143) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5375983) q[3];
sx q[3];
rz(-1.0799468) q[3];
sx q[3];
rz(-0.50118485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0665325) q[2];
sx q[2];
rz(-1.5542474) q[2];
sx q[2];
rz(1.4934322) q[2];
rz(-3.1028124) q[3];
sx q[3];
rz(-1.6438234) q[3];
sx q[3];
rz(1.405976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99603242) q[0];
sx q[0];
rz(-1.6479489) q[0];
sx q[0];
rz(2.5160568) q[0];
rz(-1.372555) q[1];
sx q[1];
rz(-0.99634975) q[1];
sx q[1];
rz(0.58403429) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82963547) q[0];
sx q[0];
rz(-2.6845508) q[0];
sx q[0];
rz(-2.5331927) q[0];
x q[1];
rz(2.4768684) q[2];
sx q[2];
rz(-1.692284) q[2];
sx q[2];
rz(0.45035502) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4702813) q[1];
sx q[1];
rz(-1.352942) q[1];
sx q[1];
rz(0.72524986) q[1];
rz(-1.8409506) q[3];
sx q[3];
rz(-0.90357257) q[3];
sx q[3];
rz(0.11303489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0293616) q[2];
sx q[2];
rz(-1.6718622) q[2];
sx q[2];
rz(1.3545797) q[2];
rz(2.3692047) q[3];
sx q[3];
rz(-1.9410917) q[3];
sx q[3];
rz(1.6105917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1188485) q[0];
sx q[0];
rz(-2.4130555) q[0];
sx q[0];
rz(1.5424315) q[0];
rz(-2.6842885) q[1];
sx q[1];
rz(-0.87304579) q[1];
sx q[1];
rz(2.9541448) q[1];
rz(1.9654033) q[2];
sx q[2];
rz(-2.9388469) q[2];
sx q[2];
rz(0.083045372) q[2];
rz(0.32021602) q[3];
sx q[3];
rz(-1.1315941) q[3];
sx q[3];
rz(-1.4109926) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
