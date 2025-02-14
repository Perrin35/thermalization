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
rz(1.4229245) q[0];
sx q[0];
rz(-2.0473502) q[0];
sx q[0];
rz(-0.25804582) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(-1.300783) q[1];
sx q[1];
rz(-2.8859477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5244183) q[0];
sx q[0];
rz(-0.99992311) q[0];
sx q[0];
rz(-2.2883795) q[0];
rz(-pi) q[1];
rz(2.7529703) q[2];
sx q[2];
rz(-1.5124694) q[2];
sx q[2];
rz(1.7494534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8430184) q[1];
sx q[1];
rz(-2.9461423) q[1];
sx q[1];
rz(-1.9324383) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37918886) q[3];
sx q[3];
rz(-2.6122836) q[3];
sx q[3];
rz(0.74931739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.51097441) q[2];
sx q[2];
rz(-1.7642085) q[2];
sx q[2];
rz(1.1154491) q[2];
rz(0.014178064) q[3];
sx q[3];
rz(-1.3445798) q[3];
sx q[3];
rz(-2.7052963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(-1.9344591) q[0];
sx q[0];
rz(-2.5196228) q[0];
sx q[0];
rz(1.2472664) q[0];
rz(-0.44218749) q[1];
sx q[1];
rz(-1.3860044) q[1];
sx q[1];
rz(-2.3242548) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38430518) q[0];
sx q[0];
rz(-1.7624859) q[0];
sx q[0];
rz(1.0377645) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27702443) q[2];
sx q[2];
rz(-1.9961341) q[2];
sx q[2];
rz(2.9556731) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4358959) q[1];
sx q[1];
rz(-1.4053939) q[1];
sx q[1];
rz(0.64343217) q[1];
rz(-pi) q[2];
rz(0.67561356) q[3];
sx q[3];
rz(-0.88425501) q[3];
sx q[3];
rz(3.0385449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3680129) q[2];
sx q[2];
rz(-0.37675884) q[2];
sx q[2];
rz(0.22029857) q[2];
rz(-1.1822654) q[3];
sx q[3];
rz(-1.9672829) q[3];
sx q[3];
rz(-0.063145414) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050506266) q[0];
sx q[0];
rz(-1.8244705) q[0];
sx q[0];
rz(1.1385981) q[0];
rz(0.41172045) q[1];
sx q[1];
rz(-2.3717334) q[1];
sx q[1];
rz(1.0134816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6665812) q[0];
sx q[0];
rz(-1.555976) q[0];
sx q[0];
rz(-2.8885452) q[0];
x q[1];
rz(-3.0092952) q[2];
sx q[2];
rz(-1.1331285) q[2];
sx q[2];
rz(-2.3486111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93597465) q[1];
sx q[1];
rz(-1.1406608) q[1];
sx q[1];
rz(2.3479152) q[1];
rz(3.0115836) q[3];
sx q[3];
rz(-0.84835669) q[3];
sx q[3];
rz(-0.60628451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98075214) q[2];
sx q[2];
rz(-1.9858805) q[2];
sx q[2];
rz(-1.8864924) q[2];
rz(0.27215019) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(-3.0009771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.960152) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(0.61035672) q[0];
rz(2.4329674) q[1];
sx q[1];
rz(-2.1253864) q[1];
sx q[1];
rz(-1.5708539) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33922903) q[0];
sx q[0];
rz(-1.7238364) q[0];
sx q[0];
rz(-1.4314851) q[0];
x q[1];
rz(0.28684692) q[2];
sx q[2];
rz(-1.45668) q[2];
sx q[2];
rz(-2.508004) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29250568) q[1];
sx q[1];
rz(-1.5226814) q[1];
sx q[1];
rz(0.27166697) q[1];
rz(0.80067486) q[3];
sx q[3];
rz(-1.281637) q[3];
sx q[3];
rz(1.2277993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2107971) q[2];
sx q[2];
rz(-1.0127298) q[2];
sx q[2];
rz(2.5861758) q[2];
rz(2.4173315) q[3];
sx q[3];
rz(-1.9925502) q[3];
sx q[3];
rz(-0.65653062) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50437462) q[0];
sx q[0];
rz(-1.1906304) q[0];
sx q[0];
rz(0.36886886) q[0];
rz(2.8952307) q[1];
sx q[1];
rz(-1.3280222) q[1];
sx q[1];
rz(-1.4341644) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.526699) q[0];
sx q[0];
rz(-1.5656359) q[0];
sx q[0];
rz(-2.3847488) q[0];
rz(1.3952012) q[2];
sx q[2];
rz(-1.1313952) q[2];
sx q[2];
rz(-2.0134157) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1885438) q[1];
sx q[1];
rz(-1.5665788) q[1];
sx q[1];
rz(-0.40178816) q[1];
x q[2];
rz(-0.09999545) q[3];
sx q[3];
rz(-1.7971562) q[3];
sx q[3];
rz(-0.55734998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4904334) q[2];
sx q[2];
rz(-1.3110524) q[2];
sx q[2];
rz(2.0595713) q[2];
rz(0.79484445) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(-1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27570462) q[0];
sx q[0];
rz(-3.0401433) q[0];
sx q[0];
rz(-0.24359447) q[0];
rz(0.98681915) q[1];
sx q[1];
rz(-2.3934264) q[1];
sx q[1];
rz(-0.93596828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0662066) q[0];
sx q[0];
rz(-0.52156007) q[0];
sx q[0];
rz(-0.88514502) q[0];
x q[1];
rz(-0.34910874) q[2];
sx q[2];
rz(-0.85452467) q[2];
sx q[2];
rz(-0.30308613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82560357) q[1];
sx q[1];
rz(-0.88615075) q[1];
sx q[1];
rz(2.0988093) q[1];
rz(-pi) q[2];
rz(-2.1706287) q[3];
sx q[3];
rz(-1.55313) q[3];
sx q[3];
rz(2.3223557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.662107) q[2];
sx q[2];
rz(-1.138849) q[2];
sx q[2];
rz(0.34995079) q[2];
rz(0.55772603) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(2.2902655) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4418942) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(-0.30174524) q[0];
rz(0.31983495) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(1.1357657) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2280884) q[0];
sx q[0];
rz(-0.64675179) q[0];
sx q[0];
rz(0.089368377) q[0];
rz(2.2280424) q[2];
sx q[2];
rz(-0.74591178) q[2];
sx q[2];
rz(-2.8532956) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8265243) q[1];
sx q[1];
rz(-0.60501912) q[1];
sx q[1];
rz(1.4304763) q[1];
rz(-3.0301827) q[3];
sx q[3];
rz(-2.3476331) q[3];
sx q[3];
rz(2.4536228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7319506) q[2];
sx q[2];
rz(-0.66764098) q[2];
sx q[2];
rz(-2.9988334) q[2];
rz(1.3879294) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9614354) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(-0.55602443) q[0];
rz(-2.9122638) q[1];
sx q[1];
rz(-1.8792968) q[1];
sx q[1];
rz(1.75846) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5189603) q[0];
sx q[0];
rz(-0.83370249) q[0];
sx q[0];
rz(-1.6181158) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11923125) q[2];
sx q[2];
rz(-0.46135739) q[2];
sx q[2];
rz(-0.90864321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.759093) q[1];
sx q[1];
rz(-1.9803932) q[1];
sx q[1];
rz(0.31633693) q[1];
x q[2];
rz(1.1874299) q[3];
sx q[3];
rz(-1.1463506) q[3];
sx q[3];
rz(-0.92680537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6250299) q[2];
sx q[2];
rz(-2.8690858) q[2];
sx q[2];
rz(1.2410835) q[2];
rz(0.062601335) q[3];
sx q[3];
rz(-1.7187985) q[3];
sx q[3];
rz(2.3144498) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0144219) q[0];
sx q[0];
rz(-1.7272471) q[0];
sx q[0];
rz(-0.14778368) q[0];
rz(-0.5087018) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(-0.37857372) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4157279) q[0];
sx q[0];
rz(-1.8047389) q[0];
sx q[0];
rz(2.7153003) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0580334) q[2];
sx q[2];
rz(-0.65750018) q[2];
sx q[2];
rz(0.12557827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3562153) q[1];
sx q[1];
rz(-1.2357724) q[1];
sx q[1];
rz(-0.32932333) q[1];
x q[2];
rz(1.4465973) q[3];
sx q[3];
rz(-0.81362766) q[3];
sx q[3];
rz(0.99909335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1754237) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(-1.8625205) q[2];
rz(1.4281979) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(-2.9960347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7487504) q[0];
sx q[0];
rz(-0.57806438) q[0];
sx q[0];
rz(0.064099126) q[0];
rz(-2.6129258) q[1];
sx q[1];
rz(-1.268498) q[1];
sx q[1];
rz(0.97506964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25697432) q[0];
sx q[0];
rz(-1.3064837) q[0];
sx q[0];
rz(-0.5075016) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3401572) q[2];
sx q[2];
rz(-0.74432997) q[2];
sx q[2];
rz(2.3542885) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2960423) q[1];
sx q[1];
rz(-0.49440372) q[1];
sx q[1];
rz(0.47459666) q[1];
rz(-0.75210877) q[3];
sx q[3];
rz(-0.32882133) q[3];
sx q[3];
rz(-3.0260835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9958682) q[2];
sx q[2];
rz(-2.4805562) q[2];
sx q[2];
rz(-2.5346942) q[2];
rz(0.25035826) q[3];
sx q[3];
rz(-1.7253877) q[3];
sx q[3];
rz(-0.50601602) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470227) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(0.49867123) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(0.90169883) q[2];
sx q[2];
rz(-0.82868262) q[2];
sx q[2];
rz(2.5862624) q[2];
rz(-3.0193083) q[3];
sx q[3];
rz(-1.3663843) q[3];
sx q[3];
rz(-2.1257675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
