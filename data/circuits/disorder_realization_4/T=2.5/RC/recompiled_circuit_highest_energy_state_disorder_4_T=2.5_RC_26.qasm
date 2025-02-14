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
rz(1.053831) q[0];
sx q[0];
rz(-2.5913936) q[0];
sx q[0];
rz(1.3681816) q[0];
rz(-0.47222459) q[1];
sx q[1];
rz(-3.0250186) q[1];
sx q[1];
rz(-0.20224686) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7305827) q[0];
sx q[0];
rz(-0.29313335) q[0];
sx q[0];
rz(-0.75384753) q[0];
rz(2.0933242) q[2];
sx q[2];
rz(-2.0179581) q[2];
sx q[2];
rz(-1.9856404) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.579548) q[1];
sx q[1];
rz(-1.288153) q[1];
sx q[1];
rz(0.94448994) q[1];
x q[2];
rz(1.6458166) q[3];
sx q[3];
rz(-2.6743691) q[3];
sx q[3];
rz(-2.644616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.54174417) q[2];
sx q[2];
rz(-0.84459633) q[2];
sx q[2];
rz(2.2005626) q[2];
rz(-2.8372676) q[3];
sx q[3];
rz(-0.86186886) q[3];
sx q[3];
rz(2.8743675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8180654) q[0];
sx q[0];
rz(-0.37779385) q[0];
sx q[0];
rz(2.3622808) q[0];
rz(1.7794973) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(1.13387) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72778217) q[0];
sx q[0];
rz(-1.6690133) q[0];
sx q[0];
rz(1.1350687) q[0];
rz(-pi) q[1];
rz(0.11793612) q[2];
sx q[2];
rz(-1.6291766) q[2];
sx q[2];
rz(-2.749325) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2782005) q[1];
sx q[1];
rz(-0.28732936) q[1];
sx q[1];
rz(-1.9598239) q[1];
rz(-1.0615223) q[3];
sx q[3];
rz(-1.5699374) q[3];
sx q[3];
rz(0.96405503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4506932) q[2];
sx q[2];
rz(-2.1096114) q[2];
sx q[2];
rz(-0.1965941) q[2];
rz(-0.96902668) q[3];
sx q[3];
rz(-1.5195547) q[3];
sx q[3];
rz(-1.903418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75152385) q[0];
sx q[0];
rz(-2.5876434) q[0];
sx q[0];
rz(0.73177904) q[0];
rz(0.36830184) q[1];
sx q[1];
rz(-0.75804561) q[1];
sx q[1];
rz(2.7751353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0822507) q[0];
sx q[0];
rz(-1.3069131) q[0];
sx q[0];
rz(-1.7543955) q[0];
rz(-pi) q[1];
rz(1.7973034) q[2];
sx q[2];
rz(-1.8069805) q[2];
sx q[2];
rz(3.0677441) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.428047) q[1];
sx q[1];
rz(-0.69880345) q[1];
sx q[1];
rz(-0.87884632) q[1];
rz(-pi) q[2];
rz(2.2499834) q[3];
sx q[3];
rz(-2.3162033) q[3];
sx q[3];
rz(-0.85435435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9590108) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(-0.89243531) q[2];
rz(-2.6871032) q[3];
sx q[3];
rz(-2.317704) q[3];
sx q[3];
rz(0.23196001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7682122) q[0];
sx q[0];
rz(-0.1425655) q[0];
sx q[0];
rz(-0.40661231) q[0];
rz(-2.3136102) q[1];
sx q[1];
rz(-1.7384638) q[1];
sx q[1];
rz(2.3060395) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0987739) q[0];
sx q[0];
rz(-1.5313498) q[0];
sx q[0];
rz(-1.7516319) q[0];
rz(-pi) q[1];
rz(-2.9391772) q[2];
sx q[2];
rz(-1.9107585) q[2];
sx q[2];
rz(2.2183403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6687739) q[1];
sx q[1];
rz(-1.7175279) q[1];
sx q[1];
rz(-2.7039578) q[1];
x q[2];
rz(-3.0441094) q[3];
sx q[3];
rz(-1.9543703) q[3];
sx q[3];
rz(2.1909379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.39997175) q[2];
sx q[2];
rz(-2.8394832) q[2];
sx q[2];
rz(0.92140222) q[2];
rz(1.7737927) q[3];
sx q[3];
rz(-0.96145815) q[3];
sx q[3];
rz(0.14149806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88857404) q[0];
sx q[0];
rz(-2.044401) q[0];
sx q[0];
rz(-0.15884037) q[0];
rz(1.5171492) q[1];
sx q[1];
rz(-0.30888638) q[1];
sx q[1];
rz(-1.812017) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8027199) q[0];
sx q[0];
rz(-3.1204528) q[0];
sx q[0];
rz(-0.52864023) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23933584) q[2];
sx q[2];
rz(-1.0221326) q[2];
sx q[2];
rz(-1.2557097) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.161474) q[1];
sx q[1];
rz(-1.9494281) q[1];
sx q[1];
rz(-2.9254854) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9333889) q[3];
sx q[3];
rz(-1.0131239) q[3];
sx q[3];
rz(0.69277908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5627731) q[2];
sx q[2];
rz(-2.4476738) q[2];
sx q[2];
rz(0.55207437) q[2];
rz(0.79365802) q[3];
sx q[3];
rz(-0.59760439) q[3];
sx q[3];
rz(-0.98646599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8296705) q[0];
sx q[0];
rz(-2.2770918) q[0];
sx q[0];
rz(-0.70575869) q[0];
rz(3.0041079) q[1];
sx q[1];
rz(-0.64422137) q[1];
sx q[1];
rz(-2.4286043) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6201545) q[0];
sx q[0];
rz(-1.7993986) q[0];
sx q[0];
rz(-2.7457353) q[0];
rz(-pi) q[1];
rz(0.51432864) q[2];
sx q[2];
rz(-1.9514958) q[2];
sx q[2];
rz(1.2027539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2299774) q[1];
sx q[1];
rz(-1.6676845) q[1];
sx q[1];
rz(1.1712892) q[1];
x q[2];
rz(-1.2014404) q[3];
sx q[3];
rz(-0.31285646) q[3];
sx q[3];
rz(-0.85858708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2116427) q[2];
sx q[2];
rz(-2.2723618) q[2];
sx q[2];
rz(-0.28406528) q[2];
rz(0.12187135) q[3];
sx q[3];
rz(-0.28212306) q[3];
sx q[3];
rz(1.4819283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26509869) q[0];
sx q[0];
rz(-2.5229186) q[0];
sx q[0];
rz(-0.70736831) q[0];
rz(-2.7012198) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(-1.7792938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61575261) q[0];
sx q[0];
rz(-2.7670015) q[0];
sx q[0];
rz(2.2994141) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84378924) q[2];
sx q[2];
rz(-1.8966881) q[2];
sx q[2];
rz(1.0032723) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8352812) q[1];
sx q[1];
rz(-1.0544485) q[1];
sx q[1];
rz(0.55037127) q[1];
rz(-2.6853722) q[3];
sx q[3];
rz(-1.4672605) q[3];
sx q[3];
rz(-0.52952784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1661487) q[2];
sx q[2];
rz(-0.41836172) q[2];
sx q[2];
rz(-2.5725906) q[2];
rz(-2.1042018) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(0.4666127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.55050945) q[0];
sx q[0];
rz(-0.35174462) q[0];
sx q[0];
rz(1.9367223) q[0];
rz(0.51271802) q[1];
sx q[1];
rz(-0.7690438) q[1];
sx q[1];
rz(2.9606294) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7962387) q[0];
sx q[0];
rz(-1.9441883) q[0];
sx q[0];
rz(0.42723421) q[0];
rz(0.34951194) q[2];
sx q[2];
rz(-1.4415359) q[2];
sx q[2];
rz(-1.845128) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9258283) q[1];
sx q[1];
rz(-1.7811532) q[1];
sx q[1];
rz(-2.1390361) q[1];
rz(-pi) q[2];
rz(3.0309903) q[3];
sx q[3];
rz(-1.6440653) q[3];
sx q[3];
rz(-1.6694809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60244954) q[2];
sx q[2];
rz(-2.4931543) q[2];
sx q[2];
rz(-3.0446206) q[2];
rz(0.55475956) q[3];
sx q[3];
rz(-1.3223038) q[3];
sx q[3];
rz(-2.7831912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6612514) q[0];
sx q[0];
rz(-2.3619409) q[0];
sx q[0];
rz(-3.0338147) q[0];
rz(-0.55094552) q[1];
sx q[1];
rz(-0.44083732) q[1];
sx q[1];
rz(-1.5239747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.658889) q[0];
sx q[0];
rz(-1.3388757) q[0];
sx q[0];
rz(-3.0085682) q[0];
rz(-pi) q[1];
rz(0.63568398) q[2];
sx q[2];
rz(-2.3909466) q[2];
sx q[2];
rz(-0.25668609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.35573449) q[1];
sx q[1];
rz(-1.5275035) q[1];
sx q[1];
rz(2.0329352) q[1];
rz(-1.6233367) q[3];
sx q[3];
rz(-1.0797622) q[3];
sx q[3];
rz(-0.72140104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23065755) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(-1.1304193) q[2];
rz(2.9039827) q[3];
sx q[3];
rz(-2.7311324) q[3];
sx q[3];
rz(2.3835278) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1199353) q[0];
sx q[0];
rz(-2.9627242) q[0];
sx q[0];
rz(0.94104952) q[0];
rz(-0.36048105) q[1];
sx q[1];
rz(-1.7345813) q[1];
sx q[1];
rz(-1.2233268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5191803) q[0];
sx q[0];
rz(-0.61953467) q[0];
sx q[0];
rz(0.34468083) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98127301) q[2];
sx q[2];
rz(-0.86893493) q[2];
sx q[2];
rz(-1.6035994) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66050038) q[1];
sx q[1];
rz(-1.5005439) q[1];
sx q[1];
rz(-2.7226647) q[1];
rz(-pi) q[2];
rz(1.2480601) q[3];
sx q[3];
rz(-0.56005037) q[3];
sx q[3];
rz(-1.7843475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2337522) q[2];
sx q[2];
rz(-1.3774104) q[2];
sx q[2];
rz(-3.0506548) q[2];
rz(-0.20445538) q[3];
sx q[3];
rz(-0.8046059) q[3];
sx q[3];
rz(-1.0819775) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37888708) q[0];
sx q[0];
rz(-1.611041) q[0];
sx q[0];
rz(1.8086717) q[0];
rz(1.7145722) q[1];
sx q[1];
rz(-2.6418229) q[1];
sx q[1];
rz(-1.5508834) q[1];
rz(-1.0574404) q[2];
sx q[2];
rz(-2.3427137) q[2];
sx q[2];
rz(-2.8165934) q[2];
rz(1.5132202) q[3];
sx q[3];
rz(-1.3575469) q[3];
sx q[3];
rz(-1.224106) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
