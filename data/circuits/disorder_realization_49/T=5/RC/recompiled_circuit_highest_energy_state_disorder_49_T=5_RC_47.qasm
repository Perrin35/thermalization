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
rz(0.81575552) q[0];
rz(-1.8176796) q[1];
sx q[1];
rz(-1.7426999) q[1];
sx q[1];
rz(0.89371347) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55339556) q[0];
sx q[0];
rz(-1.1247824) q[0];
sx q[0];
rz(1.7883401) q[0];
rz(-2.7213411) q[2];
sx q[2];
rz(-2.9300692) q[2];
sx q[2];
rz(-2.6826721) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3624139) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9245236) q[3];
sx q[3];
rz(-0.24532977) q[3];
sx q[3];
rz(-0.40553482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90014234) q[2];
sx q[2];
rz(-2.9505079) q[2];
sx q[2];
rz(-3.0639263) q[2];
rz(-0.11134722) q[3];
sx q[3];
rz(-2.1249378) q[3];
sx q[3];
rz(-2.1587423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0351008) q[0];
sx q[0];
rz(-0.74420559) q[0];
sx q[0];
rz(-2.8700854) q[0];
rz(0.60596451) q[1];
sx q[1];
rz(-0.4393591) q[1];
sx q[1];
rz(-2.4966168) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52909652) q[0];
sx q[0];
rz(-1.8865276) q[0];
sx q[0];
rz(-1.015618) q[0];
x q[1];
rz(-1.4015876) q[2];
sx q[2];
rz(-2.2300917) q[2];
sx q[2];
rz(2.7906757) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.074133679) q[1];
sx q[1];
rz(-1.3488042) q[1];
sx q[1];
rz(-1.4253474) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9927587) q[3];
sx q[3];
rz(-1.3899125) q[3];
sx q[3];
rz(2.862084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6509387) q[2];
sx q[2];
rz(-1.9927315) q[2];
sx q[2];
rz(0.33565721) q[2];
rz(3.007174) q[3];
sx q[3];
rz(-2.4498144) q[3];
sx q[3];
rz(-0.604983) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26948872) q[0];
sx q[0];
rz(-1.1815) q[0];
sx q[0];
rz(0.14542018) q[0];
rz(0.91436404) q[1];
sx q[1];
rz(-1.4375571) q[1];
sx q[1];
rz(-2.5252555) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6610044) q[0];
sx q[0];
rz(-1.0618774) q[0];
sx q[0];
rz(-0.53553217) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6583865) q[2];
sx q[2];
rz(-0.51327217) q[2];
sx q[2];
rz(-0.22944726) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65253768) q[1];
sx q[1];
rz(-1.9229865) q[1];
sx q[1];
rz(1.89025) q[1];
rz(-pi) q[2];
rz(-1.8026722) q[3];
sx q[3];
rz(-2.1926741) q[3];
sx q[3];
rz(0.14280126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99438897) q[2];
sx q[2];
rz(-0.92999593) q[2];
sx q[2];
rz(-0.27929107) q[2];
rz(0.3768557) q[3];
sx q[3];
rz(-0.91747228) q[3];
sx q[3];
rz(2.3963624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.299861) q[0];
sx q[0];
rz(-1.0693411) q[0];
sx q[0];
rz(-1.3413062) q[0];
rz(-0.74814859) q[1];
sx q[1];
rz(-2.61863) q[1];
sx q[1];
rz(-1.062692) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13650249) q[0];
sx q[0];
rz(-2.1341289) q[0];
sx q[0];
rz(-3.0045749) q[0];
rz(0.80047582) q[2];
sx q[2];
rz(-1.8006341) q[2];
sx q[2];
rz(2.2198913) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9781294) q[1];
sx q[1];
rz(-2.1390341) q[1];
sx q[1];
rz(0.37074691) q[1];
rz(-pi) q[2];
rz(-2.213988) q[3];
sx q[3];
rz(-0.91107145) q[3];
sx q[3];
rz(0.56216678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5144389) q[2];
sx q[2];
rz(-1.2320765) q[2];
sx q[2];
rz(2.3717234) q[2];
rz(1.436796) q[3];
sx q[3];
rz(-1.0154279) q[3];
sx q[3];
rz(-2.3427486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1771667) q[0];
sx q[0];
rz(-1.3097958) q[0];
sx q[0];
rz(0.31846309) q[0];
rz(-2.154649) q[1];
sx q[1];
rz(-1.8014329) q[1];
sx q[1];
rz(0.14652227) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8520235) q[0];
sx q[0];
rz(-0.98678723) q[0];
sx q[0];
rz(2.7147074) q[0];
rz(0.96168965) q[2];
sx q[2];
rz(-0.95476788) q[2];
sx q[2];
rz(0.94102695) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95899907) q[1];
sx q[1];
rz(-1.3412017) q[1];
sx q[1];
rz(0.3480546) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39590368) q[3];
sx q[3];
rz(-0.77953458) q[3];
sx q[3];
rz(0.59880873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5557308) q[2];
sx q[2];
rz(-0.18438688) q[2];
sx q[2];
rz(-1.8905852) q[2];
rz(0.74293724) q[3];
sx q[3];
rz(-1.6299959) q[3];
sx q[3];
rz(1.2062629) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10155216) q[0];
sx q[0];
rz(-1.0572301) q[0];
sx q[0];
rz(-1.2054766) q[0];
rz(2.7936753) q[1];
sx q[1];
rz(-1.9708865) q[1];
sx q[1];
rz(-1.233137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9677572) q[0];
sx q[0];
rz(-1.8223154) q[0];
sx q[0];
rz(-1.6582635) q[0];
rz(-pi) q[1];
rz(1.0516722) q[2];
sx q[2];
rz(-2.3035617) q[2];
sx q[2];
rz(-0.45853927) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72224381) q[1];
sx q[1];
rz(-2.1170685) q[1];
sx q[1];
rz(-2.6180079) q[1];
rz(-2.5830501) q[3];
sx q[3];
rz(-0.59584801) q[3];
sx q[3];
rz(2.5003025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39522383) q[2];
sx q[2];
rz(-2.2082128) q[2];
sx q[2];
rz(0.4371492) q[2];
rz(2.7002667) q[3];
sx q[3];
rz(-0.91259846) q[3];
sx q[3];
rz(0.92379409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8298892) q[0];
sx q[0];
rz(-1.403911) q[0];
sx q[0];
rz(0.29790685) q[0];
rz(0.20826134) q[1];
sx q[1];
rz(-2.2580937) q[1];
sx q[1];
rz(-2.7347402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1778616) q[0];
sx q[0];
rz(-3.047982) q[0];
sx q[0];
rz(1.9075745) q[0];
rz(-pi) q[1];
rz(-1.2804119) q[2];
sx q[2];
rz(-0.70984364) q[2];
sx q[2];
rz(0.042482201) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23644762) q[1];
sx q[1];
rz(-1.3133272) q[1];
sx q[1];
rz(-0.83579589) q[1];
x q[2];
rz(-0.86866697) q[3];
sx q[3];
rz(-2.446081) q[3];
sx q[3];
rz(1.8785541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3901861) q[2];
sx q[2];
rz(-0.69746709) q[2];
sx q[2];
rz(-0.81525272) q[2];
rz(0.32663545) q[3];
sx q[3];
rz(-1.9950461) q[3];
sx q[3];
rz(-0.013464125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1470452) q[0];
sx q[0];
rz(-2.4319686) q[0];
sx q[0];
rz(0.56353322) q[0];
rz(1.4348449) q[1];
sx q[1];
rz(-2.1610465) q[1];
sx q[1];
rz(2.837406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6582869) q[0];
sx q[0];
rz(-2.493447) q[0];
sx q[0];
rz(0.056552088) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90952342) q[2];
sx q[2];
rz(-2.3171388) q[2];
sx q[2];
rz(-2.4168991) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9619325) q[1];
sx q[1];
rz(-1.0294401) q[1];
sx q[1];
rz(2.8983745) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1564141) q[3];
sx q[3];
rz(-1.7777006) q[3];
sx q[3];
rz(-3.0690388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7131416) q[2];
sx q[2];
rz(-0.17927543) q[2];
sx q[2];
rz(0.98503867) q[2];
rz(-1.8200412) q[3];
sx q[3];
rz(-1.9165087) q[3];
sx q[3];
rz(-2.9086746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4317076) q[0];
sx q[0];
rz(-2.6490477) q[0];
sx q[0];
rz(0.4044958) q[0];
rz(-2.3996023) q[1];
sx q[1];
rz(-1.9335577) q[1];
sx q[1];
rz(-2.580339) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7405734) q[0];
sx q[0];
rz(-2.3971746) q[0];
sx q[0];
rz(0.94835001) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.77131) q[2];
sx q[2];
rz(-2.0136626) q[2];
sx q[2];
rz(2.7376769) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3041463) q[1];
sx q[1];
rz(-1.3383902) q[1];
sx q[1];
rz(-1.9327608) q[1];
x q[2];
rz(2.7453305) q[3];
sx q[3];
rz(-2.3498457) q[3];
sx q[3];
rz(1.8566673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3227417) q[2];
sx q[2];
rz(-0.35917869) q[2];
sx q[2];
rz(-0.90510577) q[2];
rz(0.94308606) q[3];
sx q[3];
rz(-1.9102996) q[3];
sx q[3];
rz(-2.4020307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1484225) q[0];
sx q[0];
rz(-0.44023308) q[0];
sx q[0];
rz(2.7711208) q[0];
rz(-0.91520339) q[1];
sx q[1];
rz(-1.7135432) q[1];
sx q[1];
rz(0.76348335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5452427) q[0];
sx q[0];
rz(-0.93726369) q[0];
sx q[0];
rz(2.0526396) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56388091) q[2];
sx q[2];
rz(-1.3016961) q[2];
sx q[2];
rz(1.0118124) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.96012989) q[1];
sx q[1];
rz(-0.69488159) q[1];
sx q[1];
rz(-1.2139266) q[1];
rz(-pi) q[2];
rz(-1.6053469) q[3];
sx q[3];
rz(-1.6240544) q[3];
sx q[3];
rz(-2.0617503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8107599) q[2];
sx q[2];
rz(-1.5786889) q[2];
sx q[2];
rz(1.2609437) q[2];
rz(1.2817945) q[3];
sx q[3];
rz(-2.1112879) q[3];
sx q[3];
rz(-1.1873881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8967165) q[0];
sx q[0];
rz(-1.5302932) q[0];
sx q[0];
rz(1.2320919) q[0];
rz(-2.2121519) q[1];
sx q[1];
rz(-0.96439958) q[1];
sx q[1];
rz(0.31698116) q[1];
rz(3.0841903) q[2];
sx q[2];
rz(-1.897503) q[2];
sx q[2];
rz(-0.31384604) q[2];
rz(0.55616641) q[3];
sx q[3];
rz(-2.5096165) q[3];
sx q[3];
rz(-1.2059042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
