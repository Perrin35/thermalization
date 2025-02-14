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
rz(-2.6925447) q[0];
sx q[0];
rz(-1.1392051) q[0];
sx q[0];
rz(2.6989302) q[0];
rz(0.17065419) q[1];
sx q[1];
rz(-2.3499188) q[1];
sx q[1];
rz(0.72905529) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013403809) q[0];
sx q[0];
rz(-0.81522757) q[0];
sx q[0];
rz(-0.56437738) q[0];
x q[1];
rz(0.38972028) q[2];
sx q[2];
rz(-1.5716296) q[2];
sx q[2];
rz(2.4454947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76293445) q[1];
sx q[1];
rz(-2.4097917) q[1];
sx q[1];
rz(-1.9373489) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7336842) q[3];
sx q[3];
rz(-2.1589557) q[3];
sx q[3];
rz(-1.5666636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.86058229) q[2];
sx q[2];
rz(-0.082110114) q[2];
sx q[2];
rz(0.78544593) q[2];
rz(-0.9663409) q[3];
sx q[3];
rz(-1.0680501) q[3];
sx q[3];
rz(0.6399703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53168374) q[0];
sx q[0];
rz(-0.56782472) q[0];
sx q[0];
rz(2.5058643) q[0];
rz(0.6334148) q[1];
sx q[1];
rz(-2.3737291) q[1];
sx q[1];
rz(-2.8107218) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71459717) q[0];
sx q[0];
rz(-1.2408371) q[0];
sx q[0];
rz(2.500588) q[0];
rz(-pi) q[1];
rz(1.517968) q[2];
sx q[2];
rz(-0.46762782) q[2];
sx q[2];
rz(-0.29948452) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2229022) q[1];
sx q[1];
rz(-1.6044334) q[1];
sx q[1];
rz(-1.6292389) q[1];
rz(1.5832354) q[3];
sx q[3];
rz(-1.6553989) q[3];
sx q[3];
rz(-1.0983262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6364381) q[2];
sx q[2];
rz(-0.81710368) q[2];
sx q[2];
rz(-0.46112296) q[2];
rz(2.4165706) q[3];
sx q[3];
rz(-2.6889763) q[3];
sx q[3];
rz(-0.66864526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.4918936) q[0];
sx q[0];
rz(-1.7576341) q[0];
sx q[0];
rz(0.53832501) q[0];
rz(0.0093983924) q[1];
sx q[1];
rz(-0.63176578) q[1];
sx q[1];
rz(-1.1993154) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0908576) q[0];
sx q[0];
rz(-1.6155302) q[0];
sx q[0];
rz(1.0629774) q[0];
rz(-0.17260562) q[2];
sx q[2];
rz(-2.4510018) q[2];
sx q[2];
rz(-1.1081072) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44502434) q[1];
sx q[1];
rz(-0.68587947) q[1];
sx q[1];
rz(0.98937757) q[1];
rz(2.8878146) q[3];
sx q[3];
rz(-2.2814085) q[3];
sx q[3];
rz(0.43076483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7678396) q[2];
sx q[2];
rz(-1.491863) q[2];
sx q[2];
rz(2.4191432) q[2];
rz(2.5433698) q[3];
sx q[3];
rz(-3.1066419) q[3];
sx q[3];
rz(-1.2178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5825321) q[0];
sx q[0];
rz(-1.8836972) q[0];
sx q[0];
rz(-3.1225358) q[0];
rz(1.4590774) q[1];
sx q[1];
rz(-0.2225114) q[1];
sx q[1];
rz(-0.51024514) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8718349) q[0];
sx q[0];
rz(-1.8352011) q[0];
sx q[0];
rz(2.1931936) q[0];
x q[1];
rz(-1.8809659) q[2];
sx q[2];
rz(-2.2686743) q[2];
sx q[2];
rz(3.0681899) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45106798) q[1];
sx q[1];
rz(-0.96123403) q[1];
sx q[1];
rz(-1.068348) q[1];
rz(-pi) q[2];
rz(0.41401569) q[3];
sx q[3];
rz(-1.7018205) q[3];
sx q[3];
rz(-1.4690635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9973008) q[2];
sx q[2];
rz(-1.4520626) q[2];
sx q[2];
rz(-2.9193381) q[2];
rz(-0.079515919) q[3];
sx q[3];
rz(-2.5955718) q[3];
sx q[3];
rz(2.5047746) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52614373) q[0];
sx q[0];
rz(-1.0960824) q[0];
sx q[0];
rz(1.3413606) q[0];
rz(-0.51097956) q[1];
sx q[1];
rz(-1.363441) q[1];
sx q[1];
rz(0.43295369) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2596672) q[0];
sx q[0];
rz(-2.8779463) q[0];
sx q[0];
rz(0.22048283) q[0];
x q[1];
rz(-1.1673981) q[2];
sx q[2];
rz(-1.7678542) q[2];
sx q[2];
rz(2.6464406) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49952415) q[1];
sx q[1];
rz(-0.41710258) q[1];
sx q[1];
rz(-2.2184371) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3687021) q[3];
sx q[3];
rz(-1.3671994) q[3];
sx q[3];
rz(-2.1264358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0114835) q[2];
sx q[2];
rz(-0.92769647) q[2];
sx q[2];
rz(0.9978655) q[2];
rz(-2.4326371) q[3];
sx q[3];
rz(-0.38781375) q[3];
sx q[3];
rz(2.1628105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3157432) q[0];
sx q[0];
rz(-2.5542927) q[0];
sx q[0];
rz(-0.89383268) q[0];
rz(-1.029344) q[1];
sx q[1];
rz(-0.7411595) q[1];
sx q[1];
rz(0.11480521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5328767) q[0];
sx q[0];
rz(-0.32793697) q[0];
sx q[0];
rz(0.84285835) q[0];
rz(-pi) q[1];
rz(2.6774339) q[2];
sx q[2];
rz(-2.4460829) q[2];
sx q[2];
rz(2.0328731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8821174) q[1];
sx q[1];
rz(-0.39813706) q[1];
sx q[1];
rz(-2.3601301) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28088681) q[3];
sx q[3];
rz(-1.6320737) q[3];
sx q[3];
rz(-2.3063502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7041695) q[2];
sx q[2];
rz(-0.88714522) q[2];
sx q[2];
rz(0.21370055) q[2];
rz(-0.63654381) q[3];
sx q[3];
rz(-0.53037363) q[3];
sx q[3];
rz(0.73591939) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8146424) q[0];
sx q[0];
rz(-2.3982168) q[0];
sx q[0];
rz(-3.0315234) q[0];
rz(-0.07668177) q[1];
sx q[1];
rz(-0.86014599) q[1];
sx q[1];
rz(-1.1383879) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0698893) q[0];
sx q[0];
rz(-1.2217297) q[0];
sx q[0];
rz(-1.3368692) q[0];
x q[1];
rz(1.5670256) q[2];
sx q[2];
rz(-2.3585883) q[2];
sx q[2];
rz(-1.4418999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0323022) q[1];
sx q[1];
rz(-1.3869973) q[1];
sx q[1];
rz(-1.079783) q[1];
x q[2];
rz(-2.8951169) q[3];
sx q[3];
rz(-1.7785871) q[3];
sx q[3];
rz(2.7026619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8818605) q[2];
sx q[2];
rz(-1.0819819) q[2];
sx q[2];
rz(-0.38121769) q[2];
rz(3.0242053) q[3];
sx q[3];
rz(-0.50506794) q[3];
sx q[3];
rz(3.05486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36670244) q[0];
sx q[0];
rz(-0.77228868) q[0];
sx q[0];
rz(0.49852398) q[0];
rz(0.82930928) q[1];
sx q[1];
rz(-2.3259951) q[1];
sx q[1];
rz(-3.0438429) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6158218) q[0];
sx q[0];
rz(-0.54057594) q[0];
sx q[0];
rz(2.1945688) q[0];
x q[1];
rz(-2.8016053) q[2];
sx q[2];
rz(-1.4302711) q[2];
sx q[2];
rz(0.7896151) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5796639) q[1];
sx q[1];
rz(-1.0439928) q[1];
sx q[1];
rz(0.71814037) q[1];
rz(-1.8923733) q[3];
sx q[3];
rz(-1.5619366) q[3];
sx q[3];
rz(1.5230501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8349614) q[2];
sx q[2];
rz(-1.9913048) q[2];
sx q[2];
rz(-0.59070307) q[2];
rz(-3.0513884) q[3];
sx q[3];
rz(-2.7019751) q[3];
sx q[3];
rz(2.2265767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-3.0331994) q[0];
sx q[0];
rz(-0.84302253) q[0];
sx q[0];
rz(-0.75476187) q[0];
rz(0.33590487) q[1];
sx q[1];
rz(-2.5867808) q[1];
sx q[1];
rz(1.0364484) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51506804) q[0];
sx q[0];
rz(-1.6881516) q[0];
sx q[0];
rz(2.4548454) q[0];
x q[1];
rz(-0.48123863) q[2];
sx q[2];
rz(-1.5889152) q[2];
sx q[2];
rz(-0.5258403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7215342) q[1];
sx q[1];
rz(-2.0871401) q[1];
sx q[1];
rz(1.5260076) q[1];
rz(-pi) q[2];
rz(-1.3410946) q[3];
sx q[3];
rz(-1.5575081) q[3];
sx q[3];
rz(2.5984327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0466517) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(-0.45647344) q[2];
rz(-0.60232317) q[3];
sx q[3];
rz(-0.18776247) q[3];
sx q[3];
rz(0.93835866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25191054) q[0];
sx q[0];
rz(-1.849181) q[0];
sx q[0];
rz(2.67814) q[0];
rz(-0.015333029) q[1];
sx q[1];
rz(-3.0749815) q[1];
sx q[1];
rz(-0.32036805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8972337) q[0];
sx q[0];
rz(-0.16288745) q[0];
sx q[0];
rz(-2.6475859) q[0];
x q[1];
rz(1.8202844) q[2];
sx q[2];
rz(-1.2659327) q[2];
sx q[2];
rz(-0.21017212) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5825504) q[1];
sx q[1];
rz(-2.1782785) q[1];
sx q[1];
rz(0.59848018) q[1];
rz(2.4418751) q[3];
sx q[3];
rz(-2.7072002) q[3];
sx q[3];
rz(0.35067973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.028712332) q[2];
sx q[2];
rz(-2.4148603) q[2];
sx q[2];
rz(-2.1325364) q[2];
rz(0.79695898) q[3];
sx q[3];
rz(-1.8877441) q[3];
sx q[3];
rz(0.33826452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0924031) q[0];
sx q[0];
rz(-1.6272463) q[0];
sx q[0];
rz(1.8820681) q[0];
rz(3.0567723) q[1];
sx q[1];
rz(-1.4823352) q[1];
sx q[1];
rz(-1.7099554) q[1];
rz(0.063719393) q[2];
sx q[2];
rz(-1.3294217) q[2];
sx q[2];
rz(0.64634993) q[2];
rz(0.30257135) q[3];
sx q[3];
rz(-1.745001) q[3];
sx q[3];
rz(0.59636084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
