OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.183455646038055) q[0];
sx q[0];
rz(2.39752105076844) q[0];
sx q[0];
rz(10.4767133951108) q[0];
rz(0.193688586354256) q[1];
sx q[1];
rz(2.50845781167085) q[1];
sx q[1];
rz(8.85177443026706) q[1];
cx q[1],q[0];
rz(-1.73180365562439) q[0];
sx q[0];
rz(2.5136725624376) q[0];
sx q[0];
rz(9.53078228085443) q[0];
rz(-2.25635147094727) q[2];
sx q[2];
rz(1.91260448296601) q[2];
sx q[2];
rz(11.9841780424039) q[2];
cx q[2],q[1];
rz(-0.515427887439728) q[1];
sx q[1];
rz(4.21607223351533) q[1];
sx q[1];
rz(12.486465907089) q[1];
rz(4.83509349822998) q[3];
sx q[3];
rz(0.700224550562449) q[3];
sx q[3];
rz(8.90007100104495) q[3];
cx q[3],q[2];
rz(-2.85945916175842) q[2];
sx q[2];
rz(2.83184749086434) q[2];
sx q[2];
rz(13.6226749181668) q[2];
rz(-3.1193962097168) q[3];
sx q[3];
rz(2.37894168694551) q[3];
sx q[3];
rz(8.0047737121503) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.231715500354767) q[0];
sx q[0];
rz(0.481193693476268) q[0];
sx q[0];
rz(9.85492241977855) q[0];
rz(-3.26928973197937) q[1];
sx q[1];
rz(5.12733558018739) q[1];
sx q[1];
rz(7.98721561431094) q[1];
cx q[1],q[0];
rz(-0.244570448994637) q[0];
sx q[0];
rz(1.87202266057069) q[0];
sx q[0];
rz(9.69488934277698) q[0];
rz(2.41434526443481) q[2];
sx q[2];
rz(3.06507160713012) q[2];
sx q[2];
rz(9.67807344197437) q[2];
cx q[2],q[1];
rz(2.39006567001343) q[1];
sx q[1];
rz(5.04030755360658) q[1];
sx q[1];
rz(8.56448290347263) q[1];
rz(-1.03112161159515) q[3];
sx q[3];
rz(-1.75753482977813) q[3];
sx q[3];
rz(13.1488461255948) q[3];
cx q[3],q[2];
rz(3.02518558502197) q[2];
sx q[2];
rz(4.6029213984781) q[2];
sx q[2];
rz(5.83775947093173) q[2];
rz(2.73235011100769) q[3];
sx q[3];
rz(5.18410411675508) q[3];
sx q[3];
rz(7.16263482569858) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.552607834339142) q[0];
sx q[0];
rz(6.19020024140412) q[0];
sx q[0];
rz(10.1395783781926) q[0];
rz(-5.22590923309326) q[1];
sx q[1];
rz(2.72092679341371) q[1];
sx q[1];
rz(12.8936333417813) q[1];
cx q[1],q[0];
rz(1.78755581378937) q[0];
sx q[0];
rz(7.68732658227021) q[0];
sx q[0];
rz(9.97338858842059) q[0];
rz(-1.13829588890076) q[2];
sx q[2];
rz(7.32786336739595) q[2];
sx q[2];
rz(7.60304305552646) q[2];
cx q[2],q[1];
rz(3.94153070449829) q[1];
sx q[1];
rz(6.13028040726716) q[1];
sx q[1];
rz(9.28909256159469) q[1];
rz(-1.82503151893616) q[3];
sx q[3];
rz(-1.2336858193106) q[3];
sx q[3];
rz(8.28084728717014) q[3];
cx q[3],q[2];
rz(-1.00322341918945) q[2];
sx q[2];
rz(5.30970898469026) q[2];
sx q[2];
rz(10.9287690877835) q[2];
rz(1.2184294462204) q[3];
sx q[3];
rz(3.5572359581762) q[3];
sx q[3];
rz(13.548381304733) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.07294356822968) q[0];
sx q[0];
rz(1.89538410504396) q[0];
sx q[0];
rz(9.26880798339053) q[0];
rz(-0.348636955022812) q[1];
sx q[1];
rz(5.67923107941682) q[1];
sx q[1];
rz(8.19177839755222) q[1];
cx q[1],q[0];
rz(-3.15916895866394) q[0];
sx q[0];
rz(4.83334711392457) q[0];
sx q[0];
rz(9.68682030438586) q[0];
rz(-1.35938596725464) q[2];
sx q[2];
rz(0.933589132624217) q[2];
sx q[2];
rz(7.87388250826999) q[2];
cx q[2],q[1];
rz(-3.02910423278809) q[1];
sx q[1];
rz(4.43073371251161) q[1];
sx q[1];
rz(6.93304679392978) q[1];
rz(1.8222668170929) q[3];
sx q[3];
rz(4.72663143475587) q[3];
sx q[3];
rz(10.2171047687451) q[3];
cx q[3],q[2];
rz(1.69845807552338) q[2];
sx q[2];
rz(-0.0360186974233905) q[2];
sx q[2];
rz(4.77173898219272) q[2];
rz(2.00214195251465) q[3];
sx q[3];
rz(4.9515359719568) q[3];
sx q[3];
rz(11.5760013818662) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.77694201469421) q[0];
sx q[0];
rz(5.86537638505036) q[0];
sx q[0];
rz(11.413348889343) q[0];
rz(-0.543683767318726) q[1];
sx q[1];
rz(5.01654973824555) q[1];
sx q[1];
rz(9.68662524818584) q[1];
cx q[1],q[0];
rz(5.51512241363525) q[0];
sx q[0];
rz(-0.0900011936849872) q[0];
sx q[0];
rz(9.80700240134403) q[0];
rz(-2.55802726745605) q[2];
sx q[2];
rz(-1.85297426382964) q[2];
sx q[2];
rz(2.50991771220371) q[2];
cx q[2],q[1];
rz(4.46022987365723) q[1];
sx q[1];
rz(4.93480363686616) q[1];
sx q[1];
rz(10.3974939942281) q[1];
rz(-2.35296893119812) q[3];
sx q[3];
rz(10.8498570044809) q[3];
sx q[3];
rz(8.71007857321903) q[3];
cx q[3],q[2];
rz(-2.52209448814392) q[2];
sx q[2];
rz(5.67749157746369) q[2];
sx q[2];
rz(7.68846461772128) q[2];
rz(-0.745534718036652) q[3];
sx q[3];
rz(8.15898052056367) q[3];
sx q[3];
rz(11.6230706930081) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.00105068227276206) q[0];
sx q[0];
rz(3.25900184561545) q[0];
sx q[0];
rz(9.77659068106815) q[0];
rz(5.84287261962891) q[1];
sx q[1];
rz(4.50706437428529) q[1];
sx q[1];
rz(6.45450208186313) q[1];
cx q[1],q[0];
rz(-2.28359436988831) q[0];
sx q[0];
rz(5.56792417367036) q[0];
sx q[0];
rz(11.3712441682737) q[0];
rz(2.32400369644165) q[2];
sx q[2];
rz(4.48701527913148) q[2];
sx q[2];
rz(4.77893016337558) q[2];
cx q[2],q[1];
rz(-0.762076675891876) q[1];
sx q[1];
rz(4.01399573882157) q[1];
sx q[1];
rz(9.39752521216079) q[1];
rz(3.83273148536682) q[3];
sx q[3];
rz(7.26558223565156) q[3];
sx q[3];
rz(12.5706581830899) q[3];
cx q[3],q[2];
rz(-6.34746932983398) q[2];
sx q[2];
rz(1.76943865616853) q[2];
sx q[2];
rz(11.6245710611264) q[2];
rz(4.29654741287231) q[3];
sx q[3];
rz(0.208081634836741) q[3];
sx q[3];
rz(11.2258546113889) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.445238590240479) q[0];
sx q[0];
rz(4.24710896809632) q[0];
sx q[0];
rz(6.28808400630161) q[0];
rz(-3.58822226524353) q[1];
sx q[1];
rz(2.54786703188951) q[1];
sx q[1];
rz(13.2543408632199) q[1];
cx q[1],q[0];
rz(-0.600583553314209) q[0];
sx q[0];
rz(4.00714418490464) q[0];
sx q[0];
rz(12.339684700958) q[0];
rz(-2.79624223709106) q[2];
sx q[2];
rz(2.16116687853868) q[2];
sx q[2];
rz(7.54699120520755) q[2];
cx q[2],q[1];
rz(-1.5060977935791) q[1];
sx q[1];
rz(0.487595470743724) q[1];
sx q[1];
rz(6.75127313136264) q[1];
rz(2.4450204372406) q[3];
sx q[3];
rz(2.13165739377076) q[3];
sx q[3];
rz(10.5323109388272) q[3];
cx q[3],q[2];
rz(-2.53043103218079) q[2];
sx q[2];
rz(0.408034237223216) q[2];
sx q[2];
rz(10.3547649145047) q[2];
rz(1.92005777359009) q[3];
sx q[3];
rz(5.17016425927217) q[3];
sx q[3];
rz(13.1597409009854) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.681948363780975) q[0];
sx q[0];
rz(4.96339562733705) q[0];
sx q[0];
rz(6.37329218386813) q[0];
rz(1.71687614917755) q[1];
sx q[1];
rz(3.7813961823755) q[1];
sx q[1];
rz(9.85986921786472) q[1];
cx q[1],q[0];
rz(2.16627073287964) q[0];
sx q[0];
rz(4.85800519784028) q[0];
sx q[0];
rz(11.4994179963987) q[0];
rz(2.23878788948059) q[2];
sx q[2];
rz(4.85495498974855) q[2];
sx q[2];
rz(6.16863820552036) q[2];
cx q[2],q[1];
rz(-0.961740911006927) q[1];
sx q[1];
rz(4.84081044991548) q[1];
sx q[1];
rz(8.80601266621753) q[1];
rz(3.14246845245361) q[3];
sx q[3];
rz(1.14643505414064) q[3];
sx q[3];
rz(13.3279294729154) q[3];
cx q[3],q[2];
rz(1.47897505760193) q[2];
sx q[2];
rz(4.2065955718332) q[2];
sx q[2];
rz(10.0509845375936) q[2];
rz(2.94467568397522) q[3];
sx q[3];
rz(5.29246154625947) q[3];
sx q[3];
rz(8.03623459338352) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.581867456436157) q[0];
sx q[0];
rz(4.83207920392091) q[0];
sx q[0];
rz(9.37049735932752) q[0];
rz(-3.56458234786987) q[1];
sx q[1];
rz(1.76616791089112) q[1];
sx q[1];
rz(11.2312005519788) q[1];
cx q[1],q[0];
rz(2.25352191925049) q[0];
sx q[0];
rz(1.0253063758188) q[0];
sx q[0];
rz(13.2130942106168) q[0];
rz(-1.21388924121857) q[2];
sx q[2];
rz(4.18165245850617) q[2];
sx q[2];
rz(11.8662726640622) q[2];
cx q[2],q[1];
rz(0.216533422470093) q[1];
sx q[1];
rz(4.19979110558564) q[1];
sx q[1];
rz(13.5260042905728) q[1];
rz(-0.459828794002533) q[3];
sx q[3];
rz(3.29526609380776) q[3];
sx q[3];
rz(10.1749544501226) q[3];
cx q[3],q[2];
rz(2.61730313301086) q[2];
sx q[2];
rz(4.60978654225404) q[2];
sx q[2];
rz(5.93168280123874) q[2];
rz(1.37378036975861) q[3];
sx q[3];
rz(3.65512994130189) q[3];
sx q[3];
rz(9.15536019801303) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.75186252593994) q[0];
sx q[0];
rz(3.4023841043287) q[0];
sx q[0];
rz(8.66561511754199) q[0];
rz(2.05793857574463) q[1];
sx q[1];
rz(7.81389823754365) q[1];
sx q[1];
rz(7.35093042849704) q[1];
cx q[1],q[0];
rz(-2.73683738708496) q[0];
sx q[0];
rz(-0.999462930364064) q[0];
sx q[0];
rz(11.5670499563138) q[0];
rz(0.277412533760071) q[2];
sx q[2];
rz(1.85458627541596) q[2];
sx q[2];
rz(14.6407861471097) q[2];
cx q[2],q[1];
rz(2.06156492233276) q[1];
sx q[1];
rz(-2.11277088324492) q[1];
sx q[1];
rz(7.37801668643161) q[1];
rz(-0.39254429936409) q[3];
sx q[3];
rz(7.01945510705049) q[3];
sx q[3];
rz(8.88797018527194) q[3];
cx q[3],q[2];
rz(1.71196734905243) q[2];
sx q[2];
rz(1.93166139920289) q[2];
sx q[2];
rz(12.7766406297605) q[2];
rz(2.35390400886536) q[3];
sx q[3];
rz(4.90039733250672) q[3];
sx q[3];
rz(9.95497403144046) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(6.31681776046753) q[0];
sx q[0];
rz(4.19154826005036) q[0];
sx q[0];
rz(8.06683025359317) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-5.56399583816528) q[1];
sx q[1];
rz(2.66473007400567) q[1];
sx q[1];
rz(14.993018603317) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-2.23522329330444) q[2];
sx q[2];
rz(2.81322443683679) q[2];
sx q[2];
rz(7.64704964160129) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-3.03346419334412) q[3];
sx q[3];
rz(1.82837644417817) q[3];
sx q[3];
rz(17.1003889799039) q[3];
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
