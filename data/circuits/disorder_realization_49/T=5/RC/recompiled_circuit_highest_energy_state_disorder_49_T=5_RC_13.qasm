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
rz(-0.10997009) q[0];
sx q[0];
rz(-2.2697544) q[0];
sx q[0];
rz(-2.3258371) q[0];
rz(-1.8176796) q[1];
sx q[1];
rz(4.5404854) q[1];
sx q[1];
rz(10.318491) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5881971) q[0];
sx q[0];
rz(-1.1247824) q[0];
sx q[0];
rz(-1.3532525) q[0];
rz(0.19359525) q[2];
sx q[2];
rz(-1.6565588) q[2];
sx q[2];
rz(2.4416358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47706931) q[1];
sx q[1];
rz(-1.1392731) q[1];
sx q[1];
rz(1.7791788) q[1];
rz(2.9245236) q[3];
sx q[3];
rz(-0.24532977) q[3];
sx q[3];
rz(0.40553482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90014234) q[2];
sx q[2];
rz(-2.9505079) q[2];
sx q[2];
rz(-3.0639263) q[2];
rz(0.11134722) q[3];
sx q[3];
rz(-2.1249378) q[3];
sx q[3];
rz(-0.98285037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0351008) q[0];
sx q[0];
rz(-2.3973871) q[0];
sx q[0];
rz(-0.27150723) q[0];
rz(0.60596451) q[1];
sx q[1];
rz(-0.4393591) q[1];
sx q[1];
rz(-2.4966168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9096268) q[0];
sx q[0];
rz(-1.0459959) q[0];
sx q[0];
rz(2.7746137) q[0];
x q[1];
rz(0.21397353) q[2];
sx q[2];
rz(-0.67751086) q[2];
sx q[2];
rz(3.0626631) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66119598) q[1];
sx q[1];
rz(-0.26473897) q[1];
sx q[1];
rz(2.5707695) q[1];
rz(-pi) q[2];
rz(-2.2521001) q[3];
sx q[3];
rz(-2.9078662) q[3];
sx q[3];
rz(0.97433486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.490654) q[2];
sx q[2];
rz(-1.9927315) q[2];
sx q[2];
rz(0.33565721) q[2];
rz(0.13441864) q[3];
sx q[3];
rz(-0.69177827) q[3];
sx q[3];
rz(2.5366096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26948872) q[0];
sx q[0];
rz(-1.9600927) q[0];
sx q[0];
rz(0.14542018) q[0];
rz(0.91436404) q[1];
sx q[1];
rz(-1.7040355) q[1];
sx q[1];
rz(2.5252555) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5443104) q[0];
sx q[0];
rz(-2.4203886) q[0];
sx q[0];
rz(-0.82996675) q[0];
rz(-0.049268289) q[2];
sx q[2];
rz(-1.0596837) q[2];
sx q[2];
rz(0.12898239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7246252) q[1];
sx q[1];
rz(-0.47096241) q[1];
sx q[1];
rz(2.4344786) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5068153) q[3];
sx q[3];
rz(-1.758681) q[3];
sx q[3];
rz(1.5646936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1472037) q[2];
sx q[2];
rz(-0.92999593) q[2];
sx q[2];
rz(0.27929107) q[2];
rz(-0.3768557) q[3];
sx q[3];
rz(-0.91747228) q[3];
sx q[3];
rz(0.74523029) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84173161) q[0];
sx q[0];
rz(-2.0722516) q[0];
sx q[0];
rz(1.8002864) q[0];
rz(-0.74814859) q[1];
sx q[1];
rz(-0.52296269) q[1];
sx q[1];
rz(1.062692) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11618075) q[0];
sx q[0];
rz(-2.5635911) q[0];
sx q[0];
rz(-1.3578174) q[0];
rz(-pi) q[1];
rz(2.8264489) q[2];
sx q[2];
rz(-0.82568554) q[2];
sx q[2];
rz(-2.7100503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1634633) q[1];
sx q[1];
rz(-1.0025585) q[1];
sx q[1];
rz(0.37074691) q[1];
rz(-pi) q[2];
x q[2];
rz(2.483401) q[3];
sx q[3];
rz(-0.88630967) q[3];
sx q[3];
rz(1.6940862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5144389) q[2];
sx q[2];
rz(-1.9095162) q[2];
sx q[2];
rz(0.76986924) q[2];
rz(-1.436796) q[3];
sx q[3];
rz(-2.1261647) q[3];
sx q[3];
rz(-2.3427486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1771667) q[0];
sx q[0];
rz(-1.3097958) q[0];
sx q[0];
rz(-2.8231296) q[0];
rz(0.98694363) q[1];
sx q[1];
rz(-1.3401597) q[1];
sx q[1];
rz(2.9950704) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6146381) q[0];
sx q[0];
rz(-1.2181158) q[0];
sx q[0];
rz(0.94278625) q[0];
rz(-pi) q[1];
rz(-2.4618838) q[2];
sx q[2];
rz(-0.83739892) q[2];
sx q[2];
rz(-1.3215166) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1720141) q[1];
sx q[1];
rz(-0.41436985) q[1];
sx q[1];
rz(2.5408391) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9349443) q[3];
sx q[3];
rz(-0.86508646) q[3];
sx q[3];
rz(-0.067506703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58586183) q[2];
sx q[2];
rz(-2.9572058) q[2];
sx q[2];
rz(-1.8905852) q[2];
rz(-0.74293724) q[3];
sx q[3];
rz(-1.6299959) q[3];
sx q[3];
rz(-1.2062629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0400405) q[0];
sx q[0];
rz(-2.0843625) q[0];
sx q[0];
rz(1.936116) q[0];
rz(2.7936753) q[1];
sx q[1];
rz(-1.9708865) q[1];
sx q[1];
rz(1.9084557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83506993) q[0];
sx q[0];
rz(-2.8756034) q[0];
sx q[0];
rz(-2.8138922) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3382931) q[2];
sx q[2];
rz(-1.9484919) q[2];
sx q[2];
rz(-0.74716728) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.55716276) q[1];
sx q[1];
rz(-2.0122156) q[1];
sx q[1];
rz(-2.1828888) q[1];
rz(-pi) q[2];
rz(0.52184409) q[3];
sx q[3];
rz(-1.2688132) q[3];
sx q[3];
rz(-2.6894231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39522383) q[2];
sx q[2];
rz(-0.93337983) q[2];
sx q[2];
rz(-2.7044435) q[2];
rz(-0.44132597) q[3];
sx q[3];
rz(-0.91259846) q[3];
sx q[3];
rz(-2.2177986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8298892) q[0];
sx q[0];
rz(-1.403911) q[0];
sx q[0];
rz(-0.29790685) q[0];
rz(0.20826134) q[1];
sx q[1];
rz(-2.2580937) q[1];
sx q[1];
rz(-2.7347402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7283467) q[0];
sx q[0];
rz(-1.6016895) q[0];
sx q[0];
rz(1.4824162) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24123206) q[2];
sx q[2];
rz(-0.89639427) q[2];
sx q[2];
rz(2.8087552) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.033483) q[1];
sx q[1];
rz(-0.86520608) q[1];
sx q[1];
rz(0.34107855) q[1];
x q[2];
rz(2.2729257) q[3];
sx q[3];
rz(-2.446081) q[3];
sx q[3];
rz(1.8785541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3901861) q[2];
sx q[2];
rz(-2.4441256) q[2];
sx q[2];
rz(0.81525272) q[2];
rz(-2.8149572) q[3];
sx q[3];
rz(-1.9950461) q[3];
sx q[3];
rz(3.1281285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1470452) q[0];
sx q[0];
rz(-2.4319686) q[0];
sx q[0];
rz(-0.56353322) q[0];
rz(-1.4348449) q[1];
sx q[1];
rz(-0.98054612) q[1];
sx q[1];
rz(2.837406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4124111) q[0];
sx q[0];
rz(-2.2177319) q[0];
sx q[0];
rz(1.6135733) q[0];
rz(-2.2772592) q[2];
sx q[2];
rz(-1.1030518) q[2];
sx q[2];
rz(-0.36004972) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9619325) q[1];
sx q[1];
rz(-1.0294401) q[1];
sx q[1];
rz(-2.8983745) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9851786) q[3];
sx q[3];
rz(-1.7777006) q[3];
sx q[3];
rz(0.072553886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7131416) q[2];
sx q[2];
rz(-0.17927543) q[2];
sx q[2];
rz(0.98503867) q[2];
rz(1.3215514) q[3];
sx q[3];
rz(-1.9165087) q[3];
sx q[3];
rz(0.2329181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4317076) q[0];
sx q[0];
rz(-2.6490477) q[0];
sx q[0];
rz(-0.4044958) q[0];
rz(2.3996023) q[1];
sx q[1];
rz(-1.9335577) q[1];
sx q[1];
rz(2.580339) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.457446) q[0];
sx q[0];
rz(-1.976891) q[0];
sx q[0];
rz(0.92828625) q[0];
rz(-2.77131) q[2];
sx q[2];
rz(-1.12793) q[2];
sx q[2];
rz(-2.7376769) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4952325) q[1];
sx q[1];
rz(-1.9226002) q[1];
sx q[1];
rz(-2.8937156) q[1];
rz(-1.1981549) q[3];
sx q[3];
rz(-0.85470457) q[3];
sx q[3];
rz(-1.3195697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81885091) q[2];
sx q[2];
rz(-2.782414) q[2];
sx q[2];
rz(-0.90510577) q[2];
rz(2.1985066) q[3];
sx q[3];
rz(-1.9102996) q[3];
sx q[3];
rz(-0.73956195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99317011) q[0];
sx q[0];
rz(-2.7013596) q[0];
sx q[0];
rz(0.37047186) q[0];
rz(-0.91520339) q[1];
sx q[1];
rz(-1.7135432) q[1];
sx q[1];
rz(-2.3781093) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5452427) q[0];
sx q[0];
rz(-0.93726369) q[0];
sx q[0];
rz(-2.0526396) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8862091) q[2];
sx q[2];
rz(-1.0295145) q[2];
sx q[2];
rz(-2.4160421) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.96012989) q[1];
sx q[1];
rz(-2.4467111) q[1];
sx q[1];
rz(1.927666) q[1];
rz(-1.6053469) q[3];
sx q[3];
rz(-1.5175382) q[3];
sx q[3];
rz(-1.0798423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3308328) q[2];
sx q[2];
rz(-1.5786889) q[2];
sx q[2];
rz(1.2609437) q[2];
rz(-1.2817945) q[3];
sx q[3];
rz(-1.0303048) q[3];
sx q[3];
rz(1.9542046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24487615) q[0];
sx q[0];
rz(-1.6112994) q[0];
sx q[0];
rz(-1.9095008) q[0];
rz(-2.2121519) q[1];
sx q[1];
rz(-0.96439958) q[1];
sx q[1];
rz(0.31698116) q[1];
rz(1.4030761) q[2];
sx q[2];
rz(-0.3315331) q[2];
sx q[2];
rz(-0.49102993) q[2];
rz(2.5852974) q[3];
sx q[3];
rz(-1.2536336) q[3];
sx q[3];
rz(-2.3118034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
