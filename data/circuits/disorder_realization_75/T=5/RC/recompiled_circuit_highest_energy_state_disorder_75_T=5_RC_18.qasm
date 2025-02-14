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
rz(-1.8259814) q[0];
sx q[0];
rz(-2.313518) q[0];
sx q[0];
rz(0.84870422) q[0];
rz(1.4915713) q[1];
sx q[1];
rz(-0.68869156) q[1];
sx q[1];
rz(0.42660776) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72683452) q[0];
sx q[0];
rz(-1.6462741) q[0];
sx q[0];
rz(-0.10087905) q[0];
x q[1];
rz(-0.83580221) q[2];
sx q[2];
rz(-2.9907132) q[2];
sx q[2];
rz(2.7250233) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6411031) q[1];
sx q[1];
rz(-1.6436542) q[1];
sx q[1];
rz(-0.44568731) q[1];
x q[2];
rz(0.60517444) q[3];
sx q[3];
rz(-0.89720336) q[3];
sx q[3];
rz(-1.0104346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2171057) q[2];
sx q[2];
rz(-1.7542398) q[2];
sx q[2];
rz(-3.1021049) q[2];
rz(-2.110179) q[3];
sx q[3];
rz(-1.02905) q[3];
sx q[3];
rz(-1.903681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51204387) q[0];
sx q[0];
rz(-1.5271674) q[0];
sx q[0];
rz(1.7426096) q[0];
rz(-1.5139187) q[1];
sx q[1];
rz(-2.2986423) q[1];
sx q[1];
rz(-1.7248076) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6607781) q[0];
sx q[0];
rz(-2.3153458) q[0];
sx q[0];
rz(-1.2673402) q[0];
rz(2.879989) q[2];
sx q[2];
rz(-0.79643476) q[2];
sx q[2];
rz(-1.4858044) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4128215) q[1];
sx q[1];
rz(-1.9660453) q[1];
sx q[1];
rz(1.4229619) q[1];
x q[2];
rz(1.6187158) q[3];
sx q[3];
rz(-1.1635492) q[3];
sx q[3];
rz(-1.4496692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.17025718) q[2];
sx q[2];
rz(-2.038326) q[2];
sx q[2];
rz(-0.7106759) q[2];
rz(-0.014287861) q[3];
sx q[3];
rz(-2.894214) q[3];
sx q[3];
rz(1.8054731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22064848) q[0];
sx q[0];
rz(-0.9592239) q[0];
sx q[0];
rz(-2.7253819) q[0];
rz(1.8572218) q[1];
sx q[1];
rz(-1.4764079) q[1];
sx q[1];
rz(-2.823337) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.153775) q[0];
sx q[0];
rz(-1.9526228) q[0];
sx q[0];
rz(-0.29819684) q[0];
x q[1];
rz(-2.2843379) q[2];
sx q[2];
rz(-0.86225715) q[2];
sx q[2];
rz(-0.015585329) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9849695) q[1];
sx q[1];
rz(-0.92069101) q[1];
sx q[1];
rz(1.7414879) q[1];
rz(0.89983799) q[3];
sx q[3];
rz(-2.7551015) q[3];
sx q[3];
rz(-0.46903601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2661813) q[2];
sx q[2];
rz(-0.77989945) q[2];
sx q[2];
rz(1.2904588) q[2];
rz(0.053704638) q[3];
sx q[3];
rz(-1.8219681) q[3];
sx q[3];
rz(-1.3857589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5829492) q[0];
sx q[0];
rz(-2.4531893) q[0];
sx q[0];
rz(2.131856) q[0];
rz(0.91521493) q[1];
sx q[1];
rz(-1.5292294) q[1];
sx q[1];
rz(-0.42542747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4998326) q[0];
sx q[0];
rz(-1.7374037) q[0];
sx q[0];
rz(-0.84979041) q[0];
rz(0.82773005) q[2];
sx q[2];
rz(-1.5166188) q[2];
sx q[2];
rz(0.82967134) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6686861) q[1];
sx q[1];
rz(-2.2232409) q[1];
sx q[1];
rz(0.5202867) q[1];
rz(-1.7336044) q[3];
sx q[3];
rz(-1.5763487) q[3];
sx q[3];
rz(2.3095363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2414744) q[2];
sx q[2];
rz(-1.590531) q[2];
sx q[2];
rz(3.0305064) q[2];
rz(2.9491718) q[3];
sx q[3];
rz(-0.40168732) q[3];
sx q[3];
rz(0.69453159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38839328) q[0];
sx q[0];
rz(-0.73047262) q[0];
sx q[0];
rz(-0.86517349) q[0];
rz(-2.6490037) q[1];
sx q[1];
rz(-1.0961696) q[1];
sx q[1];
rz(-1.7561087) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2619721) q[0];
sx q[0];
rz(-0.79085717) q[0];
sx q[0];
rz(-1.7636931) q[0];
rz(-pi) q[1];
rz(2.7215459) q[2];
sx q[2];
rz(-2.2873452) q[2];
sx q[2];
rz(1.6796215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4069479) q[1];
sx q[1];
rz(-0.67491591) q[1];
sx q[1];
rz(0.65924834) q[1];
x q[2];
rz(0.85284965) q[3];
sx q[3];
rz(-0.88969031) q[3];
sx q[3];
rz(1.8913325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8869141) q[2];
sx q[2];
rz(-2.6723537) q[2];
sx q[2];
rz(-2.4349507) q[2];
rz(0.32736579) q[3];
sx q[3];
rz(-1.2923765) q[3];
sx q[3];
rz(2.4842026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8364328) q[0];
sx q[0];
rz(-1.3494455) q[0];
sx q[0];
rz(0.35655546) q[0];
rz(0.56218475) q[1];
sx q[1];
rz(-2.0088582) q[1];
sx q[1];
rz(-0.98639375) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0537181) q[0];
sx q[0];
rz(-1.0759504) q[0];
sx q[0];
rz(-2.3943564) q[0];
rz(-0.82258921) q[2];
sx q[2];
rz(-2.6554567) q[2];
sx q[2];
rz(1.2631455) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.84556678) q[1];
sx q[1];
rz(-2.3289776) q[1];
sx q[1];
rz(-1.5350828) q[1];
rz(-pi) q[2];
rz(0.18137698) q[3];
sx q[3];
rz(-1.3311609) q[3];
sx q[3];
rz(0.58663034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76546136) q[2];
sx q[2];
rz(-0.71530801) q[2];
sx q[2];
rz(2.4410655) q[2];
rz(0.96945196) q[3];
sx q[3];
rz(-2.1078096) q[3];
sx q[3];
rz(-2.2112924) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9483865) q[0];
sx q[0];
rz(-2.8128615) q[0];
sx q[0];
rz(-2.9902003) q[0];
rz(-1.6311749) q[1];
sx q[1];
rz(-1.125234) q[1];
sx q[1];
rz(-1.6815394) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23510012) q[0];
sx q[0];
rz(-2.0564046) q[0];
sx q[0];
rz(2.4686345) q[0];
x q[1];
rz(2.7786656) q[2];
sx q[2];
rz(-2.0690326) q[2];
sx q[2];
rz(1.0403596) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5826539) q[1];
sx q[1];
rz(-2.451921) q[1];
sx q[1];
rz(-2.2037872) q[1];
x q[2];
rz(-0.93029706) q[3];
sx q[3];
rz(-0.2411763) q[3];
sx q[3];
rz(2.0335787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2363362) q[2];
sx q[2];
rz(-0.6826123) q[2];
sx q[2];
rz(0.35161099) q[2];
rz(0.44175276) q[3];
sx q[3];
rz(-0.92919246) q[3];
sx q[3];
rz(0.15171224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0996967) q[0];
sx q[0];
rz(-1.8459039) q[0];
sx q[0];
rz(-1.9805441) q[0];
rz(2.4940122) q[1];
sx q[1];
rz(-1.3787965) q[1];
sx q[1];
rz(-0.82239282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0012935972) q[0];
sx q[0];
rz(-1.9679133) q[0];
sx q[0];
rz(2.1236093) q[0];
rz(-pi) q[1];
rz(2.0089125) q[2];
sx q[2];
rz(-2.2575976) q[2];
sx q[2];
rz(-0.39114726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8257075) q[1];
sx q[1];
rz(-1.2824821) q[1];
sx q[1];
rz(2.350239) q[1];
rz(-1.137072) q[3];
sx q[3];
rz(-1.9334643) q[3];
sx q[3];
rz(0.39946242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9245727) q[2];
sx q[2];
rz(-1.9968888) q[2];
sx q[2];
rz(1.3059957) q[2];
rz(-0.71803391) q[3];
sx q[3];
rz(-0.82457232) q[3];
sx q[3];
rz(3.0927299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.7056535) q[0];
sx q[0];
rz(-2.0581364) q[0];
sx q[0];
rz(0.017024592) q[0];
rz(-1.1964993) q[1];
sx q[1];
rz(-2.7462609) q[1];
sx q[1];
rz(2.8065525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8847853) q[0];
sx q[0];
rz(-1.5809158) q[0];
sx q[0];
rz(2.1388034) q[0];
x q[1];
rz(-1.9817776) q[2];
sx q[2];
rz(-0.96945672) q[2];
sx q[2];
rz(1.3600456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6155621) q[1];
sx q[1];
rz(-2.2326734) q[1];
sx q[1];
rz(-2.9019781) q[1];
x q[2];
rz(-2.656779) q[3];
sx q[3];
rz(-1.134973) q[3];
sx q[3];
rz(-2.4143989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7177141) q[2];
sx q[2];
rz(-0.7462036) q[2];
sx q[2];
rz(2.8625873) q[2];
rz(1.772607) q[3];
sx q[3];
rz(-1.6266581) q[3];
sx q[3];
rz(2.2122808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245215) q[0];
sx q[0];
rz(-2.9947424) q[0];
sx q[0];
rz(0.76706925) q[0];
rz(2.7686139) q[1];
sx q[1];
rz(-1.6423128) q[1];
sx q[1];
rz(-2.9708718) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45324651) q[0];
sx q[0];
rz(-2.4062254) q[0];
sx q[0];
rz(0.042993025) q[0];
x q[1];
rz(2.7880048) q[2];
sx q[2];
rz(-0.43904009) q[2];
sx q[2];
rz(3.0692284) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3633299) q[1];
sx q[1];
rz(-0.40503854) q[1];
sx q[1];
rz(-0.89680775) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8324654) q[3];
sx q[3];
rz(-0.48764519) q[3];
sx q[3];
rz(1.4867099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3501222) q[2];
sx q[2];
rz(-0.56362027) q[2];
sx q[2];
rz(0.0027837022) q[2];
rz(2.2400098) q[3];
sx q[3];
rz(-1.0222579) q[3];
sx q[3];
rz(2.3367052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0881385) q[0];
sx q[0];
rz(-0.32129856) q[0];
sx q[0];
rz(-1.970485) q[0];
rz(-2.3114655) q[1];
sx q[1];
rz(-1.8142038) q[1];
sx q[1];
rz(1.289191) q[1];
rz(0.89494643) q[2];
sx q[2];
rz(-0.73113008) q[2];
sx q[2];
rz(-1.3298124) q[2];
rz(0.64408152) q[3];
sx q[3];
rz(-0.91337503) q[3];
sx q[3];
rz(2.5144284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
