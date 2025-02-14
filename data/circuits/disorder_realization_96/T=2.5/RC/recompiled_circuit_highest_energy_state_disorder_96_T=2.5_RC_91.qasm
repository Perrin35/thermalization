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
rz(1.9198298) q[0];
sx q[0];
rz(-1.3718995) q[0];
sx q[0];
rz(1.0839809) q[0];
rz(2.3594175) q[1];
sx q[1];
rz(-0.44936925) q[1];
sx q[1];
rz(2.3545797) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5840658) q[0];
sx q[0];
rz(-1.1827994) q[0];
sx q[0];
rz(0.16325133) q[0];
rz(-pi) q[1];
rz(-1.6264362) q[2];
sx q[2];
rz(-0.67650992) q[2];
sx q[2];
rz(1.8970053) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2725153) q[1];
sx q[1];
rz(-1.0306338) q[1];
sx q[1];
rz(-2.9754728) q[1];
x q[2];
rz(-0.88605864) q[3];
sx q[3];
rz(-1.3027667) q[3];
sx q[3];
rz(2.4259329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3419753) q[2];
sx q[2];
rz(-1.3222597) q[2];
sx q[2];
rz(0.70023099) q[2];
rz(1.3618943) q[3];
sx q[3];
rz(-1.2075862) q[3];
sx q[3];
rz(-0.12968682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6561061) q[0];
sx q[0];
rz(-0.34657297) q[0];
sx q[0];
rz(1.19278) q[0];
rz(-2.9876409) q[1];
sx q[1];
rz(-2.5128981) q[1];
sx q[1];
rz(-1.2492294) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5966585) q[0];
sx q[0];
rz(-1.3976263) q[0];
sx q[0];
rz(1.634266) q[0];
rz(-pi) q[1];
rz(-0.05441101) q[2];
sx q[2];
rz(-0.93428946) q[2];
sx q[2];
rz(-0.3776224) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.44012516) q[1];
sx q[1];
rz(-2.8258913) q[1];
sx q[1];
rz(-2.8786826) q[1];
rz(-1.1129679) q[3];
sx q[3];
rz(-0.71958032) q[3];
sx q[3];
rz(0.45872575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93231702) q[2];
sx q[2];
rz(-1.9220158) q[2];
sx q[2];
rz(2.2323214) q[2];
rz(-3.0026109) q[3];
sx q[3];
rz(-2.8863972) q[3];
sx q[3];
rz(2.2709258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22055498) q[0];
sx q[0];
rz(-1.9746566) q[0];
sx q[0];
rz(1.1880818) q[0];
rz(0.86337572) q[1];
sx q[1];
rz(-1.8977576) q[1];
sx q[1];
rz(-0.9695425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.506451) q[0];
sx q[0];
rz(-1.4124845) q[0];
sx q[0];
rz(-2.979605) q[0];
rz(-0.34805426) q[2];
sx q[2];
rz(-1.2173941) q[2];
sx q[2];
rz(2.5329593) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0672966) q[1];
sx q[1];
rz(-0.25416762) q[1];
sx q[1];
rz(0.39293082) q[1];
rz(-pi) q[2];
rz(-2.2754955) q[3];
sx q[3];
rz(-2.0267762) q[3];
sx q[3];
rz(2.7957819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7748176) q[2];
sx q[2];
rz(-0.25664169) q[2];
sx q[2];
rz(0.65566629) q[2];
rz(1.5244779) q[3];
sx q[3];
rz(-1.3527801) q[3];
sx q[3];
rz(0.85339439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26576385) q[0];
sx q[0];
rz(-1.1273071) q[0];
sx q[0];
rz(1.458459) q[0];
rz(-1.4153076) q[1];
sx q[1];
rz(-1.7067319) q[1];
sx q[1];
rz(-1.8728135) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44132933) q[0];
sx q[0];
rz(-2.1887795) q[0];
sx q[0];
rz(1.1412786) q[0];
rz(-0.88246302) q[2];
sx q[2];
rz(-0.75661406) q[2];
sx q[2];
rz(2.7992835) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1139196) q[1];
sx q[1];
rz(-2.259127) q[1];
sx q[1];
rz(2.9034241) q[1];
rz(-1.3425208) q[3];
sx q[3];
rz(-1.4834431) q[3];
sx q[3];
rz(-2.563208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4844674) q[2];
sx q[2];
rz(-1.7871658) q[2];
sx q[2];
rz(-1.8787059) q[2];
rz(0.2028939) q[3];
sx q[3];
rz(-2.1094567) q[3];
sx q[3];
rz(-1.3304905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48388594) q[0];
sx q[0];
rz(-1.0840451) q[0];
sx q[0];
rz(2.6935691) q[0];
rz(-1.9810642) q[1];
sx q[1];
rz(-2.0612165) q[1];
sx q[1];
rz(1.0763268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1370476) q[0];
sx q[0];
rz(-0.8595312) q[0];
sx q[0];
rz(2.8595631) q[0];
x q[1];
rz(-1.520653) q[2];
sx q[2];
rz(-1.9035982) q[2];
sx q[2];
rz(-0.34619752) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4566735) q[1];
sx q[1];
rz(-3.0701048) q[1];
sx q[1];
rz(1.0409357) q[1];
x q[2];
rz(-3.1212937) q[3];
sx q[3];
rz(-0.95385984) q[3];
sx q[3];
rz(-2.9391409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6059604) q[2];
sx q[2];
rz(-1.2781906) q[2];
sx q[2];
rz(2.7995301) q[2];
rz(2.0096807) q[3];
sx q[3];
rz(-1.071238) q[3];
sx q[3];
rz(0.38558495) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1256506) q[0];
sx q[0];
rz(-0.60369879) q[0];
sx q[0];
rz(1.1018671) q[0];
rz(-0.28533882) q[1];
sx q[1];
rz(-2.1632748) q[1];
sx q[1];
rz(2.2314821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.343276) q[0];
sx q[0];
rz(-1.7119223) q[0];
sx q[0];
rz(3.0271163) q[0];
rz(2.2533974) q[2];
sx q[2];
rz(-2.9940372) q[2];
sx q[2];
rz(0.38868586) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.460798) q[1];
sx q[1];
rz(-0.44728794) q[1];
sx q[1];
rz(1.7030925) q[1];
rz(-0.15575445) q[3];
sx q[3];
rz(-1.7475583) q[3];
sx q[3];
rz(2.8946517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3204699) q[2];
sx q[2];
rz(-1.4218825) q[2];
sx q[2];
rz(-2.087743) q[2];
rz(-2.7779135) q[3];
sx q[3];
rz(-1.6561008) q[3];
sx q[3];
rz(-0.53500879) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66154552) q[0];
sx q[0];
rz(-1.0129901) q[0];
sx q[0];
rz(0.15740982) q[0];
rz(0.62955457) q[1];
sx q[1];
rz(-1.5512356) q[1];
sx q[1];
rz(0.079364337) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.038656) q[0];
sx q[0];
rz(-3.0260822) q[0];
sx q[0];
rz(2.2282838) q[0];
x q[1];
rz(-1.7479701) q[2];
sx q[2];
rz(-1.9594155) q[2];
sx q[2];
rz(-1.7541698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2429792) q[1];
sx q[1];
rz(-0.51120394) q[1];
sx q[1];
rz(2.8752358) q[1];
rz(-pi) q[2];
rz(1.6116051) q[3];
sx q[3];
rz(-0.58678484) q[3];
sx q[3];
rz(-0.87848488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51084149) q[2];
sx q[2];
rz(-1.3151104) q[2];
sx q[2];
rz(-0.17244478) q[2];
rz(-0.6423966) q[3];
sx q[3];
rz(-0.96697092) q[3];
sx q[3];
rz(-0.35014686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.073535) q[0];
sx q[0];
rz(-1.8386766) q[0];
sx q[0];
rz(2.5750343) q[0];
rz(-0.22816518) q[1];
sx q[1];
rz(-1.4527495) q[1];
sx q[1];
rz(1.3120922) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8926548) q[0];
sx q[0];
rz(-0.83278197) q[0];
sx q[0];
rz(-2.8036462) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0242537) q[2];
sx q[2];
rz(-2.2944231) q[2];
sx q[2];
rz(0.1186419) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65441865) q[1];
sx q[1];
rz(-1.0680334) q[1];
sx q[1];
rz(0.94462402) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1817523) q[3];
sx q[3];
rz(-1.6018638) q[3];
sx q[3];
rz(-0.16060747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9606026) q[2];
sx q[2];
rz(-1.484551) q[2];
sx q[2];
rz(2.9385938) q[2];
rz(2.3024043) q[3];
sx q[3];
rz(-0.33792096) q[3];
sx q[3];
rz(0.85116974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7508271) q[0];
sx q[0];
rz(-2.6941694) q[0];
sx q[0];
rz(-0.14697337) q[0];
rz(-2.7087063) q[1];
sx q[1];
rz(-1.1914445) q[1];
sx q[1];
rz(1.8870707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0085274335) q[0];
sx q[0];
rz(-2.0957895) q[0];
sx q[0];
rz(1.9740482) q[0];
rz(1.8190838) q[2];
sx q[2];
rz(-1.5380368) q[2];
sx q[2];
rz(2.0078307) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2021823) q[1];
sx q[1];
rz(-1.9152043) q[1];
sx q[1];
rz(-0.18339524) q[1];
x q[2];
rz(0.12740429) q[3];
sx q[3];
rz(-0.16949305) q[3];
sx q[3];
rz(0.80770802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7811084) q[2];
sx q[2];
rz(-2.1238748) q[2];
sx q[2];
rz(-0.54070365) q[2];
rz(-2.8402719) q[3];
sx q[3];
rz(-0.40073985) q[3];
sx q[3];
rz(1.8436684) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28894579) q[0];
sx q[0];
rz(-1.1971373) q[0];
sx q[0];
rz(-2.4135015) q[0];
rz(2.6241265) q[1];
sx q[1];
rz(-1.6518075) q[1];
sx q[1];
rz(-0.99658406) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9835127) q[0];
sx q[0];
rz(-1.5444813) q[0];
sx q[0];
rz(1.9151494) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4696737) q[2];
sx q[2];
rz(-2.1501459) q[2];
sx q[2];
rz(-0.92321009) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4415295) q[1];
sx q[1];
rz(-2.8846488) q[1];
sx q[1];
rz(1.62613) q[1];
rz(-0.49788614) q[3];
sx q[3];
rz(-0.20789805) q[3];
sx q[3];
rz(-0.15793902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1835798) q[2];
sx q[2];
rz(-1.5637584) q[2];
sx q[2];
rz(0.070240423) q[2];
rz(3.1349414) q[3];
sx q[3];
rz(-0.17156048) q[3];
sx q[3];
rz(3.0370965) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9329279) q[0];
sx q[0];
rz(-1.7345971) q[0];
sx q[0];
rz(-1.4165184) q[0];
rz(-0.76864645) q[1];
sx q[1];
rz(-1.9192764) q[1];
sx q[1];
rz(2.8142014) q[1];
rz(-0.019975486) q[2];
sx q[2];
rz(-1.4580613) q[2];
sx q[2];
rz(-0.48403733) q[2];
rz(-0.94410789) q[3];
sx q[3];
rz(-0.44096922) q[3];
sx q[3];
rz(-0.98967688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
