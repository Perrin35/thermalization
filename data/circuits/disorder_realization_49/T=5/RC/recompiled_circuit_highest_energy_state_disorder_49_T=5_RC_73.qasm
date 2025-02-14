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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92233673) q[0];
sx q[0];
rz(-1.3748264) q[0];
sx q[0];
rz(2.6862269) q[0];
x q[1];
rz(1.4834095) q[2];
sx q[2];
rz(-1.3779216) q[2];
sx q[2];
rz(-0.88763104) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0055157) q[1];
sx q[1];
rz(-1.7598333) q[1];
sx q[1];
rz(0.4398911) q[1];
rz(-pi) q[2];
rz(-0.23979322) q[3];
sx q[3];
rz(-1.6231281) q[3];
sx q[3];
rz(0.95449846) q[3];
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
rz(0.98285037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0351008) q[0];
sx q[0];
rz(-0.74420559) q[0];
sx q[0];
rz(-0.27150723) q[0];
rz(2.5356281) q[1];
sx q[1];
rz(-2.7022336) q[1];
sx q[1];
rz(-2.4966168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52909652) q[0];
sx q[0];
rz(-1.2550651) q[0];
sx q[0];
rz(2.1259746) q[0];
rz(0.66627247) q[2];
sx q[2];
rz(-1.4372908) q[2];
sx q[2];
rz(1.324151) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.66119598) q[1];
sx q[1];
rz(-2.8768537) q[1];
sx q[1];
rz(2.5707695) q[1];
rz(-pi) q[2];
rz(-0.14883392) q[3];
sx q[3];
rz(-1.3899125) q[3];
sx q[3];
rz(-2.862084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.490654) q[2];
sx q[2];
rz(-1.9927315) q[2];
sx q[2];
rz(-2.8059354) q[2];
rz(0.13441864) q[3];
sx q[3];
rz(-2.4498144) q[3];
sx q[3];
rz(0.604983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8721039) q[0];
sx q[0];
rz(-1.1815) q[0];
sx q[0];
rz(-0.14542018) q[0];
rz(0.91436404) q[1];
sx q[1];
rz(-1.4375571) q[1];
sx q[1];
rz(0.61633715) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37164524) q[0];
sx q[0];
rz(-2.0326737) q[0];
sx q[0];
rz(-0.99528639) q[0];
x q[1];
rz(2.0824271) q[2];
sx q[2];
rz(-1.5278286) q[2];
sx q[2];
rz(1.723893) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3369477) q[1];
sx q[1];
rz(-1.8700127) q[1];
sx q[1];
rz(0.3693337) q[1];
x q[2];
rz(1.8026722) q[3];
sx q[3];
rz(-0.94891854) q[3];
sx q[3];
rz(0.14280126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99438897) q[2];
sx q[2];
rz(-0.92999593) q[2];
sx q[2];
rz(2.8623016) q[2];
rz(-2.764737) q[3];
sx q[3];
rz(-0.91747228) q[3];
sx q[3];
rz(2.3963624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.299861) q[0];
sx q[0];
rz(-1.0693411) q[0];
sx q[0];
rz(1.8002864) q[0];
rz(0.74814859) q[1];
sx q[1];
rz(-2.61863) q[1];
sx q[1];
rz(1.062692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0050902) q[0];
sx q[0];
rz(-1.0074638) q[0];
sx q[0];
rz(0.13701771) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8949365) q[2];
sx q[2];
rz(-0.79705715) q[2];
sx q[2];
rz(2.2618798) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53790435) q[1];
sx q[1];
rz(-2.474437) q[1];
sx q[1];
rz(2.0869419) q[1];
x q[2];
rz(-0.65819169) q[3];
sx q[3];
rz(-2.255283) q[3];
sx q[3];
rz(1.4475065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.62715379) q[2];
sx q[2];
rz(-1.9095162) q[2];
sx q[2];
rz(0.76986924) q[2];
rz(1.436796) q[3];
sx q[3];
rz(-1.0154279) q[3];
sx q[3];
rz(-2.3427486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1771667) q[0];
sx q[0];
rz(-1.8317969) q[0];
sx q[0];
rz(0.31846309) q[0];
rz(-0.98694363) q[1];
sx q[1];
rz(-1.3401597) q[1];
sx q[1];
rz(0.14652227) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7413899) q[0];
sx q[0];
rz(-0.70841778) q[0];
sx q[0];
rz(-1.0111459) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.429507) q[2];
sx q[2];
rz(-2.0566517) q[2];
sx q[2];
rz(0.24662247) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9695786) q[1];
sx q[1];
rz(-0.41436985) q[1];
sx q[1];
rz(0.60075356) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.745689) q[3];
sx q[3];
rz(-2.3620581) q[3];
sx q[3];
rz(-0.59880873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5557308) q[2];
sx q[2];
rz(-0.18438688) q[2];
sx q[2];
rz(-1.8905852) q[2];
rz(-0.74293724) q[3];
sx q[3];
rz(-1.5115967) q[3];
sx q[3];
rz(1.2062629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.10155216) q[0];
sx q[0];
rz(-1.0572301) q[0];
sx q[0];
rz(-1.936116) q[0];
rz(2.7936753) q[1];
sx q[1];
rz(-1.9708865) q[1];
sx q[1];
rz(1.9084557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7228112) q[0];
sx q[0];
rz(-1.486088) q[0];
sx q[0];
rz(-2.8891488) q[0];
x q[1];
rz(0.50384028) q[2];
sx q[2];
rz(-0.8692534) q[2];
sx q[2];
rz(1.1654677) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55716276) q[1];
sx q[1];
rz(-1.129377) q[1];
sx q[1];
rz(2.1828888) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55854256) q[3];
sx q[3];
rz(-2.5457446) q[3];
sx q[3];
rz(2.5003025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39522383) q[2];
sx q[2];
rz(-2.2082128) q[2];
sx q[2];
rz(0.4371492) q[2];
rz(-0.44132597) q[3];
sx q[3];
rz(-0.91259846) q[3];
sx q[3];
rz(0.92379409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8298892) q[0];
sx q[0];
rz(-1.403911) q[0];
sx q[0];
rz(0.29790685) q[0];
rz(-0.20826134) q[1];
sx q[1];
rz(-2.2580937) q[1];
sx q[1];
rz(2.7347402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16028743) q[0];
sx q[0];
rz(-1.6591342) q[0];
sx q[0];
rz(-3.1105785) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88201875) q[2];
sx q[2];
rz(-1.758496) q[2];
sx q[2];
rz(-1.7512013) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.905145) q[1];
sx q[1];
rz(-1.3133272) q[1];
sx q[1];
rz(2.3057968) q[1];
x q[2];
rz(2.6471944) q[3];
sx q[3];
rz(-2.0819849) q[3];
sx q[3];
rz(2.7123812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75140658) q[2];
sx q[2];
rz(-2.4441256) q[2];
sx q[2];
rz(0.81525272) q[2];
rz(0.32663545) q[3];
sx q[3];
rz(-1.1465466) q[3];
sx q[3];
rz(0.013464125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470452) q[0];
sx q[0];
rz(-0.70962405) q[0];
sx q[0];
rz(2.5780594) q[0];
rz(-1.7067478) q[1];
sx q[1];
rz(-2.1610465) q[1];
sx q[1];
rz(2.837406) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6582869) q[0];
sx q[0];
rz(-2.493447) q[0];
sx q[0];
rz(-0.056552088) q[0];
rz(-pi) q[1];
rz(2.5553914) q[2];
sx q[2];
rz(-2.1888141) q[2];
sx q[2];
rz(1.5780592) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1796602) q[1];
sx q[1];
rz(-2.1121526) q[1];
sx q[1];
rz(-0.24321817) q[1];
rz(1.9851786) q[3];
sx q[3];
rz(-1.7777006) q[3];
sx q[3];
rz(3.0690388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4284511) q[2];
sx q[2];
rz(-2.9623172) q[2];
sx q[2];
rz(2.156554) q[2];
rz(-1.3215514) q[3];
sx q[3];
rz(-1.225084) q[3];
sx q[3];
rz(0.2329181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4317076) q[0];
sx q[0];
rz(-0.49254492) q[0];
sx q[0];
rz(-0.4044958) q[0];
rz(0.74199039) q[1];
sx q[1];
rz(-1.208035) q[1];
sx q[1];
rz(-0.56125364) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9674751) q[0];
sx q[0];
rz(-0.9878648) q[0];
sx q[0];
rz(2.6487104) q[0];
x q[1];
rz(-2.2225631) q[2];
sx q[2];
rz(-2.5723151) q[2];
sx q[2];
rz(-2.809466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8617647) q[1];
sx q[1];
rz(-2.7142254) q[1];
sx q[1];
rz(0.98160569) q[1];
x q[2];
rz(2.3901479) q[3];
sx q[3];
rz(-1.2925694) q[3];
sx q[3];
rz(-3.1415526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3227417) q[2];
sx q[2];
rz(-0.35917869) q[2];
sx q[2];
rz(0.90510577) q[2];
rz(2.1985066) q[3];
sx q[3];
rz(-1.9102996) q[3];
sx q[3];
rz(2.4020307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(2.2263893) q[1];
sx q[1];
rz(-1.4280495) q[1];
sx q[1];
rz(2.3781093) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32578308) q[0];
sx q[0];
rz(-1.9535583) q[0];
sx q[0];
rz(0.69212691) q[0];
rz(2.5777117) q[2];
sx q[2];
rz(-1.3016961) q[2];
sx q[2];
rz(2.1297803) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88958538) q[1];
sx q[1];
rz(-1.3452054) q[1];
sx q[1];
rz(0.90771339) q[1];
rz(-pi) q[2];
rz(-3.0883028) q[3];
sx q[3];
rz(-1.6052979) q[3];
sx q[3];
rz(2.6487987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8107599) q[2];
sx q[2];
rz(-1.5786889) q[2];
sx q[2];
rz(1.8806489) q[2];
rz(1.2817945) q[3];
sx q[3];
rz(-1.0303048) q[3];
sx q[3];
rz(-1.9542046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8967165) q[0];
sx q[0];
rz(-1.6112994) q[0];
sx q[0];
rz(-1.9095008) q[0];
rz(-2.2121519) q[1];
sx q[1];
rz(-0.96439958) q[1];
sx q[1];
rz(0.31698116) q[1];
rz(1.8980044) q[2];
sx q[2];
rz(-1.6251593) q[2];
sx q[2];
rz(-1.9030824) q[2];
rz(-0.55616641) q[3];
sx q[3];
rz(-0.63197613) q[3];
sx q[3];
rz(1.9356884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
