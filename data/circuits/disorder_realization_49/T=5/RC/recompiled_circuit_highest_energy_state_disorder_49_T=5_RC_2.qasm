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
rz(0.87183824) q[0];
sx q[0];
rz(11.750615) q[0];
rz(-1.8176796) q[1];
sx q[1];
rz(-1.7426999) q[1];
sx q[1];
rz(0.89371347) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1146671) q[0];
sx q[0];
rz(-0.4930149) q[0];
sx q[0];
rz(0.42401029) q[0];
rz(-pi) q[1];
rz(-1.4834095) q[2];
sx q[2];
rz(-1.7636711) q[2];
sx q[2];
rz(-0.88763104) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1360769) q[1];
sx q[1];
rz(-1.7598333) q[1];
sx q[1];
rz(2.7017016) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2170691) q[3];
sx q[3];
rz(-2.8962629) q[3];
sx q[3];
rz(-0.40553482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2414503) q[2];
sx q[2];
rz(-0.19108471) q[2];
sx q[2];
rz(-3.0639263) q[2];
rz(3.0302454) q[3];
sx q[3];
rz(-1.0166549) q[3];
sx q[3];
rz(2.1587423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0351008) q[0];
sx q[0];
rz(-0.74420559) q[0];
sx q[0];
rz(-0.27150723) q[0];
rz(-0.60596451) q[1];
sx q[1];
rz(-2.7022336) q[1];
sx q[1];
rz(-2.4966168) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9096268) q[0];
sx q[0];
rz(-2.0955968) q[0];
sx q[0];
rz(0.36697893) q[0];
x q[1];
rz(-0.21397353) q[2];
sx q[2];
rz(-2.4640818) q[2];
sx q[2];
rz(3.0626631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4803967) q[1];
sx q[1];
rz(-2.8768537) q[1];
sx q[1];
rz(-2.5707695) q[1];
rz(-pi) q[2];
rz(-0.88949253) q[3];
sx q[3];
rz(-2.9078662) q[3];
sx q[3];
rz(-0.97433486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.490654) q[2];
sx q[2];
rz(-1.9927315) q[2];
sx q[2];
rz(2.8059354) q[2];
rz(3.007174) q[3];
sx q[3];
rz(-0.69177827) q[3];
sx q[3];
rz(-2.5366096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8721039) q[0];
sx q[0];
rz(-1.9600927) q[0];
sx q[0];
rz(-2.9961725) q[0];
rz(-0.91436404) q[1];
sx q[1];
rz(-1.4375571) q[1];
sx q[1];
rz(-0.61633715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37164524) q[0];
sx q[0];
rz(-1.108919) q[0];
sx q[0];
rz(0.99528639) q[0];
rz(-pi) q[1];
rz(1.6583865) q[2];
sx q[2];
rz(-0.51327217) q[2];
sx q[2];
rz(-2.9121454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7246252) q[1];
sx q[1];
rz(-0.47096241) q[1];
sx q[1];
rz(-0.70711406) q[1];
x q[2];
rz(-0.63477739) q[3];
sx q[3];
rz(-1.3829117) q[3];
sx q[3];
rz(-1.5646936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1472037) q[2];
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
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.299861) q[0];
sx q[0];
rz(-2.0722516) q[0];
sx q[0];
rz(-1.8002864) q[0];
rz(-0.74814859) q[1];
sx q[1];
rz(-0.52296269) q[1];
sx q[1];
rz(1.062692) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5077909) q[0];
sx q[0];
rz(-1.4550545) q[0];
sx q[0];
rz(-1.0032038) q[0];
rz(-1.8949365) q[2];
sx q[2];
rz(-2.3445355) q[2];
sx q[2];
rz(-2.2618798) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53790435) q[1];
sx q[1];
rz(-2.474437) q[1];
sx q[1];
rz(-1.0546507) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92760466) q[3];
sx q[3];
rz(-0.91107145) q[3];
sx q[3];
rz(-0.56216678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62715379) q[2];
sx q[2];
rz(-1.2320765) q[2];
sx q[2];
rz(2.3717234) q[2];
rz(-1.436796) q[3];
sx q[3];
rz(-2.1261647) q[3];
sx q[3];
rz(0.79884401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.964426) q[0];
sx q[0];
rz(-1.3097958) q[0];
sx q[0];
rz(2.8231296) q[0];
rz(0.98694363) q[1];
sx q[1];
rz(-1.3401597) q[1];
sx q[1];
rz(-0.14652227) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40020271) q[0];
sx q[0];
rz(-0.70841778) q[0];
sx q[0];
rz(2.1304467) q[0];
rz(-pi) q[1];
rz(0.96168965) q[2];
sx q[2];
rz(-2.1868248) q[2];
sx q[2];
rz(-0.94102695) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9695786) q[1];
sx q[1];
rz(-0.41436985) q[1];
sx q[1];
rz(0.60075356) q[1];
x q[2];
rz(-2.4022465) q[3];
sx q[3];
rz(-1.845318) q[3];
sx q[3];
rz(1.2609466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58586183) q[2];
sx q[2];
rz(-2.9572058) q[2];
sx q[2];
rz(-1.2510074) q[2];
rz(-0.74293724) q[3];
sx q[3];
rz(-1.5115967) q[3];
sx q[3];
rz(-1.9353297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0400405) q[0];
sx q[0];
rz(-1.0572301) q[0];
sx q[0];
rz(1.936116) q[0];
rz(-0.34791738) q[1];
sx q[1];
rz(-1.9708865) q[1];
sx q[1];
rz(1.9084557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83506993) q[0];
sx q[0];
rz(-0.26598922) q[0];
sx q[0];
rz(-0.32770043) q[0];
x q[1];
rz(-2.6377524) q[2];
sx q[2];
rz(-0.8692534) q[2];
sx q[2];
rz(1.1654677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5603139) q[1];
sx q[1];
rz(-2.403861) q[1];
sx q[1];
rz(-2.2590265) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55854256) q[3];
sx q[3];
rz(-0.59584801) q[3];
sx q[3];
rz(-2.5003025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39522383) q[2];
sx q[2];
rz(-2.2082128) q[2];
sx q[2];
rz(2.7044435) q[2];
rz(-2.7002667) q[3];
sx q[3];
rz(-2.2289942) q[3];
sx q[3];
rz(0.92379409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31170347) q[0];
sx q[0];
rz(-1.403911) q[0];
sx q[0];
rz(2.8436858) q[0];
rz(-0.20826134) q[1];
sx q[1];
rz(-0.88349897) q[1];
sx q[1];
rz(-2.7347402) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1778616) q[0];
sx q[0];
rz(-0.093610659) q[0];
sx q[0];
rz(1.9075745) q[0];
rz(-pi) q[1];
rz(-1.2804119) q[2];
sx q[2];
rz(-0.70984364) q[2];
sx q[2];
rz(0.042482201) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23644762) q[1];
sx q[1];
rz(-1.8282654) q[1];
sx q[1];
rz(-2.3057968) q[1];
rz(-1.0034542) q[3];
sx q[3];
rz(-1.1441244) q[3];
sx q[3];
rz(-0.88374352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3901861) q[2];
sx q[2];
rz(-2.4441256) q[2];
sx q[2];
rz(-2.3263399) q[2];
rz(-0.32663545) q[3];
sx q[3];
rz(-1.1465466) q[3];
sx q[3];
rz(-0.013464125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470452) q[0];
sx q[0];
rz(-0.70962405) q[0];
sx q[0];
rz(-0.56353322) q[0];
rz(1.7067478) q[1];
sx q[1];
rz(-0.98054612) q[1];
sx q[1];
rz(-0.30418667) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2741843) q[0];
sx q[0];
rz(-1.5366669) q[0];
sx q[0];
rz(-2.4942168) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5553914) q[2];
sx q[2];
rz(-0.95277856) q[2];
sx q[2];
rz(1.5780592) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9619325) q[1];
sx q[1];
rz(-1.0294401) q[1];
sx q[1];
rz(2.8983745) q[1];
rz(1.0902152) q[3];
sx q[3];
rz(-2.6811185) q[3];
sx q[3];
rz(-1.0613031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4284511) q[2];
sx q[2];
rz(-0.17927543) q[2];
sx q[2];
rz(2.156554) q[2];
rz(-1.8200412) q[3];
sx q[3];
rz(-1.9165087) q[3];
sx q[3];
rz(0.2329181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.709885) q[0];
sx q[0];
rz(-0.49254492) q[0];
sx q[0];
rz(-0.4044958) q[0];
rz(2.3996023) q[1];
sx q[1];
rz(-1.9335577) q[1];
sx q[1];
rz(2.580339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40101926) q[0];
sx q[0];
rz(-0.74441806) q[0];
sx q[0];
rz(2.1932426) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0414342) q[2];
sx q[2];
rz(-1.9038892) q[2];
sx q[2];
rz(-1.3317219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2798279) q[1];
sx q[1];
rz(-0.4273673) q[1];
sx q[1];
rz(-2.159987) q[1];
rz(2.7453305) q[3];
sx q[3];
rz(-2.3498457) q[3];
sx q[3];
rz(-1.2849253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3227417) q[2];
sx q[2];
rz(-2.782414) q[2];
sx q[2];
rz(0.90510577) q[2];
rz(0.94308606) q[3];
sx q[3];
rz(-1.9102996) q[3];
sx q[3];
rz(0.73956195) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1484225) q[0];
sx q[0];
rz(-0.44023308) q[0];
sx q[0];
rz(0.37047186) q[0];
rz(-0.91520339) q[1];
sx q[1];
rz(-1.4280495) q[1];
sx q[1];
rz(-0.76348335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8158096) q[0];
sx q[0];
rz(-1.9535583) q[0];
sx q[0];
rz(0.69212691) q[0];
rz(-0.56388091) q[2];
sx q[2];
rz(-1.8398966) q[2];
sx q[2];
rz(-2.1297803) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88958538) q[1];
sx q[1];
rz(-1.7963872) q[1];
sx q[1];
rz(0.90771339) q[1];
rz(-pi) q[2];
rz(-1.6053469) q[3];
sx q[3];
rz(-1.6240544) q[3];
sx q[3];
rz(-2.0617503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8107599) q[2];
sx q[2];
rz(-1.5629038) q[2];
sx q[2];
rz(1.8806489) q[2];
rz(1.2817945) q[3];
sx q[3];
rz(-1.0303048) q[3];
sx q[3];
rz(1.1873881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
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
rz(2.2121519) q[1];
sx q[1];
rz(-2.1771931) q[1];
sx q[1];
rz(-2.8246115) q[1];
rz(1.2435883) q[2];
sx q[2];
rz(-1.5164334) q[2];
sx q[2];
rz(1.2385102) q[2];
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
