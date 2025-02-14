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
rz(1.323913) q[1];
sx q[1];
rz(-1.3988928) q[1];
sx q[1];
rz(-0.89371347) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55339556) q[0];
sx q[0];
rz(-1.1247824) q[0];
sx q[0];
rz(1.3532525) q[0];
rz(-pi) q[1];
rz(2.9479974) q[2];
sx q[2];
rz(-1.6565588) q[2];
sx q[2];
rz(-2.4416358) q[2];
x q[3];
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
rz(-1.3624139) q[1];
rz(-pi) q[2];
rz(-1.6246666) q[3];
sx q[3];
rz(-1.3313378) q[3];
sx q[3];
rz(2.5125063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90014234) q[2];
sx q[2];
rz(-2.9505079) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10649189) q[0];
sx q[0];
rz(-2.3973871) q[0];
sx q[0];
rz(2.8700854) q[0];
rz(0.60596451) q[1];
sx q[1];
rz(-0.4393591) q[1];
sx q[1];
rz(-2.4966168) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52909652) q[0];
sx q[0];
rz(-1.2550651) q[0];
sx q[0];
rz(2.1259746) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21397353) q[2];
sx q[2];
rz(-2.4640818) q[2];
sx q[2];
rz(-3.0626631) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6771705) q[1];
sx q[1];
rz(-1.4289411) q[1];
sx q[1];
rz(-2.9173098) q[1];
x q[2];
rz(-2.9927587) q[3];
sx q[3];
rz(-1.7516802) q[3];
sx q[3];
rz(0.27950865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.490654) q[2];
sx q[2];
rz(-1.1488612) q[2];
sx q[2];
rz(-0.33565721) q[2];
rz(-3.007174) q[3];
sx q[3];
rz(-2.4498144) q[3];
sx q[3];
rz(0.604983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26948872) q[0];
sx q[0];
rz(-1.9600927) q[0];
sx q[0];
rz(0.14542018) q[0];
rz(-2.2272286) q[1];
sx q[1];
rz(-1.4375571) q[1];
sx q[1];
rz(0.61633715) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7699474) q[0];
sx q[0];
rz(-2.0326737) q[0];
sx q[0];
rz(-2.1463063) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0591655) q[2];
sx q[2];
rz(-1.613764) q[2];
sx q[2];
rz(-1.4176996) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3369477) q[1];
sx q[1];
rz(-1.27158) q[1];
sx q[1];
rz(2.772259) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3389204) q[3];
sx q[3];
rz(-0.94891854) q[3];
sx q[3];
rz(-0.14280126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1472037) q[2];
sx q[2];
rz(-0.92999593) q[2];
sx q[2];
rz(-2.8623016) q[2];
rz(-0.3768557) q[3];
sx q[3];
rz(-2.2241204) q[3];
sx q[3];
rz(-0.74523029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.299861) q[0];
sx q[0];
rz(-2.0722516) q[0];
sx q[0];
rz(1.8002864) q[0];
rz(2.3934441) q[1];
sx q[1];
rz(-0.52296269) q[1];
sx q[1];
rz(-2.0789007) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5077909) q[0];
sx q[0];
rz(-1.4550545) q[0];
sx q[0];
rz(1.0032038) q[0];
rz(1.8949365) q[2];
sx q[2];
rz(-2.3445355) q[2];
sx q[2];
rz(2.2618798) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9781294) q[1];
sx q[1];
rz(-1.0025585) q[1];
sx q[1];
rz(2.7708457) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65819169) q[3];
sx q[3];
rz(-0.88630967) q[3];
sx q[3];
rz(1.6940862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5144389) q[2];
sx q[2];
rz(-1.9095162) q[2];
sx q[2];
rz(2.3717234) q[2];
rz(1.7047966) q[3];
sx q[3];
rz(-2.1261647) q[3];
sx q[3];
rz(-2.3427486) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
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
rz(-1.8014329) q[1];
sx q[1];
rz(-0.14652227) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7413899) q[0];
sx q[0];
rz(-2.4331749) q[0];
sx q[0];
rz(2.1304467) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96168965) q[2];
sx q[2];
rz(-2.1868248) q[2];
sx q[2];
rz(2.2005657) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1825936) q[1];
sx q[1];
rz(-1.800391) q[1];
sx q[1];
rz(0.3480546) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9349443) q[3];
sx q[3];
rz(-0.86508646) q[3];
sx q[3];
rz(0.067506703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.58586183) q[2];
sx q[2];
rz(-0.18438688) q[2];
sx q[2];
rz(-1.2510074) q[2];
rz(-2.3986554) q[3];
sx q[3];
rz(-1.6299959) q[3];
sx q[3];
rz(1.2062629) q[3];
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
rz(-pi) q[0];
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
rz(-1.233137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83506993) q[0];
sx q[0];
rz(-0.26598922) q[0];
sx q[0];
rz(-2.8138922) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50384028) q[2];
sx q[2];
rz(-0.8692534) q[2];
sx q[2];
rz(-1.9761249) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5603139) q[1];
sx q[1];
rz(-0.73773161) q[1];
sx q[1];
rz(-0.88256617) q[1];
x q[2];
rz(1.9157655) q[3];
sx q[3];
rz(-1.0747834) q[3];
sx q[3];
rz(1.2880052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7463688) q[2];
sx q[2];
rz(-2.2082128) q[2];
sx q[2];
rz(2.7044435) q[2];
rz(-2.7002667) q[3];
sx q[3];
rz(-0.91259846) q[3];
sx q[3];
rz(-0.92379409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8298892) q[0];
sx q[0];
rz(-1.7376816) q[0];
sx q[0];
rz(2.8436858) q[0];
rz(-0.20826134) q[1];
sx q[1];
rz(-2.2580937) q[1];
sx q[1];
rz(-0.40685245) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9813052) q[0];
sx q[0];
rz(-1.4824585) q[0];
sx q[0];
rz(-3.1105785) q[0];
rz(-pi) q[1];
rz(0.88201875) q[2];
sx q[2];
rz(-1.758496) q[2];
sx q[2];
rz(1.3903914) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5326724) q[1];
sx q[1];
rz(-2.3708276) q[1];
sx q[1];
rz(1.9449598) q[1];
rz(2.2729257) q[3];
sx q[3];
rz(-2.446081) q[3];
sx q[3];
rz(-1.2630386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3901861) q[2];
sx q[2];
rz(-0.69746709) q[2];
sx q[2];
rz(0.81525272) q[2];
rz(0.32663545) q[3];
sx q[3];
rz(-1.9950461) q[3];
sx q[3];
rz(3.1281285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99454749) q[0];
sx q[0];
rz(-0.70962405) q[0];
sx q[0];
rz(-2.5780594) q[0];
rz(1.4348449) q[1];
sx q[1];
rz(-0.98054612) q[1];
sx q[1];
rz(0.30418667) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4124111) q[0];
sx q[0];
rz(-0.92386073) q[0];
sx q[0];
rz(1.5280194) q[0];
rz(-pi) q[1];
rz(0.86433347) q[2];
sx q[2];
rz(-1.1030518) q[2];
sx q[2];
rz(2.7815429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4106928) q[1];
sx q[1];
rz(-0.58847702) q[1];
sx q[1];
rz(-1.1898349) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0513774) q[3];
sx q[3];
rz(-2.6811185) q[3];
sx q[3];
rz(1.0613031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7131416) q[2];
sx q[2];
rz(-0.17927543) q[2];
sx q[2];
rz(-0.98503867) q[2];
rz(1.8200412) q[3];
sx q[3];
rz(-1.225084) q[3];
sx q[3];
rz(-2.9086746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.709885) q[0];
sx q[0];
rz(-0.49254492) q[0];
sx q[0];
rz(-0.4044958) q[0];
rz(-0.74199039) q[1];
sx q[1];
rz(-1.208035) q[1];
sx q[1];
rz(0.56125364) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7405734) q[0];
sx q[0];
rz(-0.74441806) q[0];
sx q[0];
rz(0.94835001) q[0];
rz(-0.37028266) q[2];
sx q[2];
rz(-1.12793) q[2];
sx q[2];
rz(-0.40391573) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.64636014) q[1];
sx q[1];
rz(-1.2189924) q[1];
sx q[1];
rz(2.8937156) q[1];
rz(2.7453305) q[3];
sx q[3];
rz(-2.3498457) q[3];
sx q[3];
rz(1.8566673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3227417) q[2];
sx q[2];
rz(-0.35917869) q[2];
sx q[2];
rz(-2.2364869) q[2];
rz(0.94308606) q[3];
sx q[3];
rz(-1.2312931) q[3];
sx q[3];
rz(-0.73956195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99317011) q[0];
sx q[0];
rz(-0.44023308) q[0];
sx q[0];
rz(0.37047186) q[0];
rz(0.91520339) q[1];
sx q[1];
rz(-1.7135432) q[1];
sx q[1];
rz(-0.76348335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32578308) q[0];
sx q[0];
rz(-1.1880344) q[0];
sx q[0];
rz(-0.69212691) q[0];
x q[1];
rz(2.6652226) q[2];
sx q[2];
rz(-2.5231276) q[2];
sx q[2];
rz(-0.16099421) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50824805) q[1];
sx q[1];
rz(-0.92735433) q[1];
sx q[1];
rz(0.28336759) q[1];
rz(-3.0883028) q[3];
sx q[3];
rz(-1.5362947) q[3];
sx q[3];
rz(-2.6487987) q[3];
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
rz(1.8597982) q[3];
sx q[3];
rz(-1.0303048) q[3];
sx q[3];
rz(-1.1873881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24487615) q[0];
sx q[0];
rz(-1.6112994) q[0];
sx q[0];
rz(-1.9095008) q[0];
rz(2.2121519) q[1];
sx q[1];
rz(-2.1771931) q[1];
sx q[1];
rz(-2.8246115) q[1];
rz(-1.7385165) q[2];
sx q[2];
rz(-0.3315331) q[2];
sx q[2];
rz(-0.49102993) q[2];
rz(1.2019602) q[3];
sx q[3];
rz(-2.0963674) q[3];
sx q[3];
rz(2.5921303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
