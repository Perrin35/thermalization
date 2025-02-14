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
rz(1.2827058) q[0];
sx q[0];
rz(-2.1252706) q[0];
sx q[0];
rz(1.8860201) q[0];
rz(-1.3610871) q[1];
sx q[1];
rz(-0.95870107) q[1];
sx q[1];
rz(0.60671848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5131305) q[0];
sx q[0];
rz(-2.2431122) q[0];
sx q[0];
rz(-0.99654128) q[0];
rz(-pi) q[1];
rz(1.8357359) q[2];
sx q[2];
rz(-1.5393672) q[2];
sx q[2];
rz(-0.72223483) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9470939) q[1];
sx q[1];
rz(-1.7003807) q[1];
sx q[1];
rz(1.2034077) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6912724) q[3];
sx q[3];
rz(-1.7306149) q[3];
sx q[3];
rz(-0.32630703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9386193) q[2];
sx q[2];
rz(-1.9846658) q[2];
sx q[2];
rz(0.17091664) q[2];
rz(1.1275229) q[3];
sx q[3];
rz(-2.5850962) q[3];
sx q[3];
rz(0.053430406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.663986) q[0];
sx q[0];
rz(-0.32653102) q[0];
sx q[0];
rz(-0.68847454) q[0];
rz(1.7636048) q[1];
sx q[1];
rz(-2.1207899) q[1];
sx q[1];
rz(-0.14150208) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48442885) q[0];
sx q[0];
rz(-1.1714897) q[0];
sx q[0];
rz(-2.5998373) q[0];
x q[1];
rz(-0.42948855) q[2];
sx q[2];
rz(-1.6613243) q[2];
sx q[2];
rz(1.2538101) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0558171) q[1];
sx q[1];
rz(-1.9250828) q[1];
sx q[1];
rz(1.8721168) q[1];
rz(-pi) q[2];
rz(-1.3745295) q[3];
sx q[3];
rz(-0.39159039) q[3];
sx q[3];
rz(2.151769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3839533) q[2];
sx q[2];
rz(-2.6944104) q[2];
sx q[2];
rz(0.028707061) q[2];
rz(0.38255295) q[3];
sx q[3];
rz(-1.5111204) q[3];
sx q[3];
rz(0.20146519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(0.89556995) q[0];
sx q[0];
rz(-2.0992794) q[0];
sx q[0];
rz(-2.7699455) q[0];
rz(0.21013513) q[1];
sx q[1];
rz(-0.34966436) q[1];
sx q[1];
rz(1.9336611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71825224) q[0];
sx q[0];
rz(-2.1499942) q[0];
sx q[0];
rz(0.96286003) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5919331) q[2];
sx q[2];
rz(-1.653228) q[2];
sx q[2];
rz(-1.7475413) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74151285) q[1];
sx q[1];
rz(-0.99435213) q[1];
sx q[1];
rz(2.4680572) q[1];
x q[2];
rz(-1.5299093) q[3];
sx q[3];
rz(-1.2256116) q[3];
sx q[3];
rz(0.57905771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5189884) q[2];
sx q[2];
rz(-3.0486139) q[2];
sx q[2];
rz(0.63817111) q[2];
rz(-2.7291164) q[3];
sx q[3];
rz(-1.4700438) q[3];
sx q[3];
rz(1.1023869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07688044) q[0];
sx q[0];
rz(-2.2314254) q[0];
sx q[0];
rz(1.8338715) q[0];
rz(-0.67963302) q[1];
sx q[1];
rz(-0.72703397) q[1];
sx q[1];
rz(-1.9576498) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9480431) q[0];
sx q[0];
rz(-1.8947161) q[0];
sx q[0];
rz(-2.4452423) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6865191) q[2];
sx q[2];
rz(-1.3365796) q[2];
sx q[2];
rz(-0.17143347) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.013742208) q[1];
sx q[1];
rz(-0.35687414) q[1];
sx q[1];
rz(2.9452717) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77149646) q[3];
sx q[3];
rz(-2.960254) q[3];
sx q[3];
rz(2.4490631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1188941) q[2];
sx q[2];
rz(-2.9939632) q[2];
sx q[2];
rz(1.0559319) q[2];
rz(2.8583156) q[3];
sx q[3];
rz(-1.2669468) q[3];
sx q[3];
rz(-1.383708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0530171) q[0];
sx q[0];
rz(-0.96240369) q[0];
sx q[0];
rz(-1.4085294) q[0];
rz(-0.22964302) q[1];
sx q[1];
rz(-0.85669986) q[1];
sx q[1];
rz(-1.1913258) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8386779) q[0];
sx q[0];
rz(-1.9698799) q[0];
sx q[0];
rz(-0.15388925) q[0];
x q[1];
rz(3.0158305) q[2];
sx q[2];
rz(-1.874141) q[2];
sx q[2];
rz(-1.6605008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3670259) q[1];
sx q[1];
rz(-1.3210249) q[1];
sx q[1];
rz(2.7109409) q[1];
rz(-2.2226589) q[3];
sx q[3];
rz(-2.0410547) q[3];
sx q[3];
rz(1.4138427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64629897) q[2];
sx q[2];
rz(-0.36094347) q[2];
sx q[2];
rz(1.6634644) q[2];
rz(-2.9804001) q[3];
sx q[3];
rz(-1.5321923) q[3];
sx q[3];
rz(-0.78981367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5021055) q[0];
sx q[0];
rz(-0.4244856) q[0];
sx q[0];
rz(0.91833997) q[0];
rz(2.8847252) q[1];
sx q[1];
rz(-0.99595064) q[1];
sx q[1];
rz(2.0253983) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74280069) q[0];
sx q[0];
rz(-2.0896119) q[0];
sx q[0];
rz(-2.1543845) q[0];
rz(-1.6565422) q[2];
sx q[2];
rz(-1.298279) q[2];
sx q[2];
rz(2.9086824) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0829385) q[1];
sx q[1];
rz(-1.8095329) q[1];
sx q[1];
rz(-1.3707364) q[1];
rz(-pi) q[2];
rz(-0.243354) q[3];
sx q[3];
rz(-1.4019308) q[3];
sx q[3];
rz(-0.7272184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9548698) q[2];
sx q[2];
rz(-2.335151) q[2];
sx q[2];
rz(2.7346129) q[2];
rz(-0.61378971) q[3];
sx q[3];
rz(-1.611462) q[3];
sx q[3];
rz(0.49155244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3864022) q[0];
sx q[0];
rz(-0.84119216) q[0];
sx q[0];
rz(0.92591539) q[0];
rz(2.6689802) q[1];
sx q[1];
rz(-2.0874529) q[1];
sx q[1];
rz(-0.47017631) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0553207) q[0];
sx q[0];
rz(-2.4013858) q[0];
sx q[0];
rz(1.6946778) q[0];
rz(2.2165259) q[2];
sx q[2];
rz(-1.888141) q[2];
sx q[2];
rz(0.44337526) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9035569) q[1];
sx q[1];
rz(-0.78552526) q[1];
sx q[1];
rz(3.0901647) q[1];
rz(-pi) q[2];
rz(-2.3005465) q[3];
sx q[3];
rz(-1.5980835) q[3];
sx q[3];
rz(2.5430665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4818695) q[2];
sx q[2];
rz(-2.0707776) q[2];
sx q[2];
rz(-0.83079633) q[2];
rz(-3.0610436) q[3];
sx q[3];
rz(-2.060067) q[3];
sx q[3];
rz(2.1848333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.047121) q[0];
sx q[0];
rz(-1.6166649) q[0];
sx q[0];
rz(-2.9834874) q[0];
rz(2.415601) q[1];
sx q[1];
rz(-0.43015614) q[1];
sx q[1];
rz(-0.37011883) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67664941) q[0];
sx q[0];
rz(-1.9029362) q[0];
sx q[0];
rz(-1.667883) q[0];
x q[1];
rz(2.3838777) q[2];
sx q[2];
rz(-1.9279459) q[2];
sx q[2];
rz(-2.0983608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8565048) q[1];
sx q[1];
rz(-2.5743432) q[1];
sx q[1];
rz(-2.3598266) q[1];
rz(-pi) q[2];
rz(-1.0256938) q[3];
sx q[3];
rz(-1.2152367) q[3];
sx q[3];
rz(-2.4073383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98081723) q[2];
sx q[2];
rz(-1.846401) q[2];
sx q[2];
rz(0.13713169) q[2];
rz(-1.6454227) q[3];
sx q[3];
rz(-2.4479595) q[3];
sx q[3];
rz(-0.48367286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5982323) q[0];
sx q[0];
rz(-1.0022663) q[0];
sx q[0];
rz(-0.12635669) q[0];
rz(-2.5670746) q[1];
sx q[1];
rz(-0.9205598) q[1];
sx q[1];
rz(-0.7472907) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4536323) q[0];
sx q[0];
rz(-1.9575141) q[0];
sx q[0];
rz(-2.9403015) q[0];
x q[1];
rz(2.4805675) q[2];
sx q[2];
rz(-2.3213561) q[2];
sx q[2];
rz(0.30414061) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14366985) q[1];
sx q[1];
rz(-1.2005245) q[1];
sx q[1];
rz(-1.0940432) q[1];
rz(-pi) q[2];
rz(0.46183698) q[3];
sx q[3];
rz(-0.29134068) q[3];
sx q[3];
rz(2.9614465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0390465) q[2];
sx q[2];
rz(-2.0206082) q[2];
sx q[2];
rz(-2.5940564) q[2];
rz(1.3012137) q[3];
sx q[3];
rz(-2.0042714) q[3];
sx q[3];
rz(3.014452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.309677) q[0];
sx q[0];
rz(-1.3818106) q[0];
sx q[0];
rz(-0.58050138) q[0];
rz(0.25513395) q[1];
sx q[1];
rz(-2.3105123) q[1];
sx q[1];
rz(-1.1855804) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0017173926) q[0];
sx q[0];
rz(-2.1809362) q[0];
sx q[0];
rz(-0.82307477) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9337148) q[2];
sx q[2];
rz(-1.5693773) q[2];
sx q[2];
rz(-0.05677536) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8557786) q[1];
sx q[1];
rz(-1.4818945) q[1];
sx q[1];
rz(-2.6707021) q[1];
rz(-0.2106726) q[3];
sx q[3];
rz(-0.86189358) q[3];
sx q[3];
rz(-1.7514995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4127976) q[2];
sx q[2];
rz(-2.2788861) q[2];
sx q[2];
rz(-2.1346788) q[2];
rz(-1.5236731) q[3];
sx q[3];
rz(-2.1606052) q[3];
sx q[3];
rz(-1.1716918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11378743) q[0];
sx q[0];
rz(-1.2600949) q[0];
sx q[0];
rz(-2.8366198) q[0];
rz(1.2420568) q[1];
sx q[1];
rz(-1.2868953) q[1];
sx q[1];
rz(-2.3931265) q[1];
rz(1.2075049) q[2];
sx q[2];
rz(-0.45798326) q[2];
sx q[2];
rz(1.9275673) q[2];
rz(2.2119568) q[3];
sx q[3];
rz(-0.84992483) q[3];
sx q[3];
rz(0.58569943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
