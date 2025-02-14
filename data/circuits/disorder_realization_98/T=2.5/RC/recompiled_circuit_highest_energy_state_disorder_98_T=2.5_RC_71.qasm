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
rz(-0.12914817) q[0];
sx q[0];
rz(-1.669786) q[0];
sx q[0];
rz(0.9653402) q[0];
rz(0.020429285) q[1];
sx q[1];
rz(-1.483622) q[1];
sx q[1];
rz(-1.3774011) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55104461) q[0];
sx q[0];
rz(-0.88860308) q[0];
sx q[0];
rz(2.5709573) q[0];
x q[1];
rz(-1.5280058) q[2];
sx q[2];
rz(-2.4952609) q[2];
sx q[2];
rz(0.94852322) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5038576) q[1];
sx q[1];
rz(-1.4252932) q[1];
sx q[1];
rz(2.9521431) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3404779) q[3];
sx q[3];
rz(-1.5571888) q[3];
sx q[3];
rz(2.1134562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91780245) q[2];
sx q[2];
rz(-1.6518355) q[2];
sx q[2];
rz(-1.6415143) q[2];
rz(-0.7817868) q[3];
sx q[3];
rz(-2.2512071) q[3];
sx q[3];
rz(-2.3894943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016945275) q[0];
sx q[0];
rz(-0.77777672) q[0];
sx q[0];
rz(-2.1710904) q[0];
rz(-1.3264725) q[1];
sx q[1];
rz(-1.7992203) q[1];
sx q[1];
rz(2.2296947) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0835595) q[0];
sx q[0];
rz(-2.5290997) q[0];
sx q[0];
rz(0.29077323) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3746337) q[2];
sx q[2];
rz(-1.7387783) q[2];
sx q[2];
rz(-2.2998435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40373793) q[1];
sx q[1];
rz(-0.7016088) q[1];
sx q[1];
rz(1.6340294) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1158783) q[3];
sx q[3];
rz(-1.5568241) q[3];
sx q[3];
rz(-0.25532237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9087002) q[2];
sx q[2];
rz(-1.1129271) q[2];
sx q[2];
rz(-0.98928893) q[2];
rz(-0.42012897) q[3];
sx q[3];
rz(-1.0003072) q[3];
sx q[3];
rz(-2.4205128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99783889) q[0];
sx q[0];
rz(-1.1622575) q[0];
sx q[0];
rz(1.0379399) q[0];
rz(-2.2833917) q[1];
sx q[1];
rz(-2.61519) q[1];
sx q[1];
rz(-0.16955489) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35992453) q[0];
sx q[0];
rz(-1.5667652) q[0];
sx q[0];
rz(1.4272593) q[0];
x q[1];
rz(1.1716258) q[2];
sx q[2];
rz(-1.8695651) q[2];
sx q[2];
rz(-2.5577161) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9537064) q[1];
sx q[1];
rz(-2.0665209) q[1];
sx q[1];
rz(0.33457054) q[1];
rz(-pi) q[2];
rz(1.9334458) q[3];
sx q[3];
rz(-1.7154105) q[3];
sx q[3];
rz(0.1000769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9618591) q[2];
sx q[2];
rz(-2.6229975) q[2];
sx q[2];
rz(2.7355984) q[2];
rz(0.77754846) q[3];
sx q[3];
rz(-2.8113139) q[3];
sx q[3];
rz(-2.3269261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6430214) q[0];
sx q[0];
rz(-1.850147) q[0];
sx q[0];
rz(2.2136069) q[0];
rz(3.0827177) q[1];
sx q[1];
rz(-1.6935655) q[1];
sx q[1];
rz(1.4307129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3979523) q[0];
sx q[0];
rz(-2.7317606) q[0];
sx q[0];
rz(1.3751956) q[0];
rz(1.769999) q[2];
sx q[2];
rz(-1.2143233) q[2];
sx q[2];
rz(1.3717331) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3763172) q[1];
sx q[1];
rz(-2.5936454) q[1];
sx q[1];
rz(2.0085232) q[1];
rz(-pi) q[2];
rz(-0.30117463) q[3];
sx q[3];
rz(-2.1780925) q[3];
sx q[3];
rz(-3.0454896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4579939) q[2];
sx q[2];
rz(-1.4983838) q[2];
sx q[2];
rz(1.9963473) q[2];
rz(-0.84732071) q[3];
sx q[3];
rz(-1.3342131) q[3];
sx q[3];
rz(1.5020717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.1640846) q[0];
sx q[0];
rz(-1.9590398) q[0];
sx q[0];
rz(-2.7567647) q[0];
rz(-0.79967868) q[1];
sx q[1];
rz(-2.622602) q[1];
sx q[1];
rz(-1.5934561) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.528397) q[0];
sx q[0];
rz(-1.3269182) q[0];
sx q[0];
rz(0.24922483) q[0];
x q[1];
rz(1.1742915) q[2];
sx q[2];
rz(-1.7616211) q[2];
sx q[2];
rz(0.79939524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5729294) q[1];
sx q[1];
rz(-1.8844114) q[1];
sx q[1];
rz(0.36797087) q[1];
rz(1.0043275) q[3];
sx q[3];
rz(-0.62823717) q[3];
sx q[3];
rz(0.35997501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3758214) q[2];
sx q[2];
rz(-2.4479726) q[2];
sx q[2];
rz(-0.080246298) q[2];
rz(-2.5806228) q[3];
sx q[3];
rz(-1.4813981) q[3];
sx q[3];
rz(-0.62275732) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65492594) q[0];
sx q[0];
rz(-2.4764562) q[0];
sx q[0];
rz(1.8810062) q[0];
rz(-0.27443019) q[1];
sx q[1];
rz(-1.471289) q[1];
sx q[1];
rz(2.4226277) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.270416) q[0];
sx q[0];
rz(-0.76379062) q[0];
sx q[0];
rz(-1.0538573) q[0];
rz(-pi) q[1];
rz(0.0098835398) q[2];
sx q[2];
rz(-0.82429574) q[2];
sx q[2];
rz(3.1205999) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5605041) q[1];
sx q[1];
rz(-1.6916654) q[1];
sx q[1];
rz(-2.493119) q[1];
rz(-2.0893065) q[3];
sx q[3];
rz(-2.720359) q[3];
sx q[3];
rz(-0.13430922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6545973) q[2];
sx q[2];
rz(-1.8751112) q[2];
sx q[2];
rz(-2.5385762) q[2];
rz(-1.4740137) q[3];
sx q[3];
rz(-2.1173729) q[3];
sx q[3];
rz(-2.730864) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5176158) q[0];
sx q[0];
rz(-1.6083953) q[0];
sx q[0];
rz(0.18727592) q[0];
rz(-2.1693443) q[1];
sx q[1];
rz(-0.15521237) q[1];
sx q[1];
rz(2.7211199) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0213326) q[0];
sx q[0];
rz(-1.9870523) q[0];
sx q[0];
rz(0.72159213) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0275317) q[2];
sx q[2];
rz(-1.8993371) q[2];
sx q[2];
rz(1.7787041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5125888) q[1];
sx q[1];
rz(-1.2930451) q[1];
sx q[1];
rz(-1.1729097) q[1];
rz(-pi) q[2];
rz(1.2240846) q[3];
sx q[3];
rz(-0.60255749) q[3];
sx q[3];
rz(1.6438705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.095852701) q[2];
sx q[2];
rz(-2.0407712) q[2];
sx q[2];
rz(-2.8908758) q[2];
rz(-1.2480674) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(2.5417476) q[3];
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
rz(pi/2) q[0];
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
rz(-1.8641758) q[0];
sx q[0];
rz(-0.59159788) q[0];
sx q[0];
rz(-0.94171062) q[0];
rz(-0.038453728) q[1];
sx q[1];
rz(-1.7105303) q[1];
sx q[1];
rz(-3.0252735) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.250688) q[0];
sx q[0];
rz(-0.90849344) q[0];
sx q[0];
rz(-2.566889) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8031249) q[2];
sx q[2];
rz(-0.69455494) q[2];
sx q[2];
rz(1.6603927) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6814057) q[1];
sx q[1];
rz(-1.7509772) q[1];
sx q[1];
rz(2.8799675) q[1];
rz(-pi) q[2];
rz(-2.7523858) q[3];
sx q[3];
rz(-1.4688244) q[3];
sx q[3];
rz(1.9736279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2945127) q[2];
sx q[2];
rz(-0.81081644) q[2];
sx q[2];
rz(2.1152367) q[2];
rz(-1.9319084) q[3];
sx q[3];
rz(-1.0299725) q[3];
sx q[3];
rz(-2.1799555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.472979) q[0];
sx q[0];
rz(-1.0740148) q[0];
sx q[0];
rz(0.37856722) q[0];
rz(1.1096795) q[1];
sx q[1];
rz(-2.8329284) q[1];
sx q[1];
rz(3.1139156) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97229703) q[0];
sx q[0];
rz(-2.1624273) q[0];
sx q[0];
rz(-1.1280498) q[0];
x q[1];
rz(-1.354486) q[2];
sx q[2];
rz(-1.3900818) q[2];
sx q[2];
rz(-1.2265151) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31071404) q[1];
sx q[1];
rz(-2.7355477) q[1];
sx q[1];
rz(-1.0881523) q[1];
x q[2];
rz(-2.1975193) q[3];
sx q[3];
rz(-2.6289796) q[3];
sx q[3];
rz(-2.1533898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3045584) q[2];
sx q[2];
rz(-2.1874032) q[2];
sx q[2];
rz(2.7216116) q[2];
rz(-0.84152451) q[3];
sx q[3];
rz(-2.1244815) q[3];
sx q[3];
rz(-2.7628472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751752) q[0];
sx q[0];
rz(-0.11688047) q[0];
sx q[0];
rz(0.28512678) q[0];
rz(2.9410703) q[1];
sx q[1];
rz(-1.9384117) q[1];
sx q[1];
rz(0.83555046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0203945) q[0];
sx q[0];
rz(-0.9011974) q[0];
sx q[0];
rz(-0.59359896) q[0];
x q[1];
rz(0.64705683) q[2];
sx q[2];
rz(-1.5655883) q[2];
sx q[2];
rz(-0.13089779) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3491464) q[1];
sx q[1];
rz(-2.0172152) q[1];
sx q[1];
rz(-1.2013776) q[1];
x q[2];
rz(0.23563175) q[3];
sx q[3];
rz(-1.6088369) q[3];
sx q[3];
rz(0.63562449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42651132) q[2];
sx q[2];
rz(-2.0501523) q[2];
sx q[2];
rz(0.69768989) q[2];
rz(-0.52608025) q[3];
sx q[3];
rz(-2.3487921) q[3];
sx q[3];
rz(1.5066159) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1220916) q[0];
sx q[0];
rz(-2.2631336) q[0];
sx q[0];
rz(2.7418131) q[0];
rz(-1.9922235) q[1];
sx q[1];
rz(-2.414357) q[1];
sx q[1];
rz(2.3946708) q[1];
rz(-2.7928945) q[2];
sx q[2];
rz(-0.92338466) q[2];
sx q[2];
rz(0.79104214) q[2];
rz(2.5935843) q[3];
sx q[3];
rz(-1.3153793) q[3];
sx q[3];
rz(1.8438189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
