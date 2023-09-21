OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(-1.707466) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2375862) q[0];
sx q[0];
rz(-1.307784) q[0];
sx q[0];
rz(-1.0884652) q[0];
x q[1];
rz(0.11169545) q[2];
sx q[2];
rz(-1.7061491) q[2];
sx q[2];
rz(0.40458194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0093706) q[1];
sx q[1];
rz(-1.8144061) q[1];
sx q[1];
rz(-1.1019215) q[1];
x q[2];
rz(-2.9207346) q[3];
sx q[3];
rz(-2.0005895) q[3];
sx q[3];
rz(-2.5151099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.819954) q[2];
sx q[2];
rz(-1.7724089) q[2];
sx q[2];
rz(2.3036172) q[2];
rz(-2.6485802) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-0.078991927) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3943966) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.863742) q[0];
rz(-0.17678075) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(0.4321672) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55789253) q[0];
sx q[0];
rz(-2.0895134) q[0];
sx q[0];
rz(-1.8899263) q[0];
rz(-pi) q[1];
rz(-0.82168174) q[2];
sx q[2];
rz(-0.5558008) q[2];
sx q[2];
rz(-2.1205714) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4914815) q[1];
sx q[1];
rz(-1.1307798) q[1];
sx q[1];
rz(0.12932175) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0057955) q[3];
sx q[3];
rz(-1.5917935) q[3];
sx q[3];
rz(2.3137623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5923578) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(0.48669997) q[2];
rz(-1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(-0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040722672) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(-1.4886645) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(-2.6352077) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9837512) q[0];
sx q[0];
rz(-0.39641532) q[0];
sx q[0];
rz(2.0435964) q[0];
rz(-pi) q[1];
x q[1];
rz(1.330553) q[2];
sx q[2];
rz(-0.85075399) q[2];
sx q[2];
rz(-2.4400997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.10889974) q[1];
sx q[1];
rz(-1.1316205) q[1];
sx q[1];
rz(-1.9119309) q[1];
x q[2];
rz(2.0154325) q[3];
sx q[3];
rz(-2.7362842) q[3];
sx q[3];
rz(2.0733548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4425519) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(0.59147269) q[2];
rz(2.5555723) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(-1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1699003) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(-2.3024094) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(1.5485839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73491053) q[0];
sx q[0];
rz(-1.4797749) q[0];
sx q[0];
rz(-2.0216366) q[0];
x q[1];
rz(0.84393878) q[2];
sx q[2];
rz(-0.91274777) q[2];
sx q[2];
rz(-2.7968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5316481) q[1];
sx q[1];
rz(-2.2713486) q[1];
sx q[1];
rz(-1.0611666) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2876835) q[3];
sx q[3];
rz(-1.1554171) q[3];
sx q[3];
rz(1.1113885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.30535355) q[2];
sx q[2];
rz(-2.2576015) q[2];
sx q[2];
rz(0.099686064) q[2];
rz(2.1827407) q[3];
sx q[3];
rz(-1.3189664) q[3];
sx q[3];
rz(1.3249741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56617671) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(-2.8856522) q[0];
rz(-2.6804965) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(-0.76006132) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81239031) q[0];
sx q[0];
rz(-0.15604067) q[0];
sx q[0];
rz(0.72326707) q[0];
rz(0.063644479) q[2];
sx q[2];
rz(-1.9846989) q[2];
sx q[2];
rz(2.6457583) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6872245) q[1];
sx q[1];
rz(-1.6974653) q[1];
sx q[1];
rz(1.441799) q[1];
x q[2];
rz(0.36565904) q[3];
sx q[3];
rz(-0.95698157) q[3];
sx q[3];
rz(-0.97096503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5710859) q[2];
sx q[2];
rz(-1.3723624) q[2];
sx q[2];
rz(2.467353) q[2];
rz(-0.21480602) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-3.1242483) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53133416) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(-1.9624788) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(-2.0746322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93766312) q[0];
sx q[0];
rz(-1.5563037) q[0];
sx q[0];
rz(0.020676215) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7166769) q[2];
sx q[2];
rz(-1.9746466) q[2];
sx q[2];
rz(2.0770819) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21642099) q[1];
sx q[1];
rz(-0.74837084) q[1];
sx q[1];
rz(-1.5009576) q[1];
x q[2];
rz(-0.79640572) q[3];
sx q[3];
rz(-1.4561597) q[3];
sx q[3];
rz(0.20533268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.431488) q[2];
sx q[2];
rz(-1.2892712) q[2];
sx q[2];
rz(0.55994326) q[2];
rz(-2.4152749) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85754919) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(-1.42111) q[0];
rz(0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(0.85817671) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35408033) q[0];
sx q[0];
rz(-1.4424099) q[0];
sx q[0];
rz(-1.7187198) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.981295) q[2];
sx q[2];
rz(-1.5294642) q[2];
sx q[2];
rz(1.2683887) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3837636) q[1];
sx q[1];
rz(-1.5997636) q[1];
sx q[1];
rz(-1.5555192) q[1];
x q[2];
rz(-2.6303597) q[3];
sx q[3];
rz(-2.5066262) q[3];
sx q[3];
rz(-0.2583897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8537366) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(1.8161592) q[2];
rz(-2.251513) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(-0.98852283) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778075) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(-0.37471399) q[0];
rz(0.97887865) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(-1.7920866) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6357248) q[0];
sx q[0];
rz(-0.74698193) q[0];
sx q[0];
rz(2.5542459) q[0];
x q[1];
rz(-1.5478412) q[2];
sx q[2];
rz(-0.85283961) q[2];
sx q[2];
rz(-0.1644451) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0950349) q[1];
sx q[1];
rz(-1.3890424) q[1];
sx q[1];
rz(-3.0912116) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8324864) q[3];
sx q[3];
rz(-1.040254) q[3];
sx q[3];
rz(1.7390651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65138856) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(-0.46869579) q[2];
rz(-1.9474585) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(1.3635925) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63012183) q[0];
sx q[0];
rz(-2.2352495) q[0];
sx q[0];
rz(-1.2619031) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(1.6040241) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9984734) q[0];
sx q[0];
rz(-2.7612918) q[0];
sx q[0];
rz(2.2024676) q[0];
rz(1.4666918) q[2];
sx q[2];
rz(-1.9874007) q[2];
sx q[2];
rz(-2.5660851) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9201587) q[1];
sx q[1];
rz(-2.4281574) q[1];
sx q[1];
rz(2.9677797) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92026199) q[3];
sx q[3];
rz(-0.55763054) q[3];
sx q[3];
rz(-2.695431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.039915446) q[2];
sx q[2];
rz(-0.95255178) q[2];
sx q[2];
rz(2.3802479) q[2];
rz(0.90041655) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.802357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(2.9246869) q[0];
rz(2.5096109) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(2.1868618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2749598) q[0];
sx q[0];
rz(-1.2506335) q[0];
sx q[0];
rz(-0.76036705) q[0];
rz(1.9913313) q[2];
sx q[2];
rz(-2.0805801) q[2];
sx q[2];
rz(1.2812986) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6001544) q[1];
sx q[1];
rz(-2.6791875) q[1];
sx q[1];
rz(-2.9701783) q[1];
x q[2];
rz(2.6984152) q[3];
sx q[3];
rz(-2.157353) q[3];
sx q[3];
rz(-3.0122258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1251936) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(2.1968502) q[2];
rz(2.7567806) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(-2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5007297) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(-2.2568933) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(-2.2434071) q[2];
sx q[2];
rz(-2.2834416) q[2];
sx q[2];
rz(-1.7590547) q[2];
rz(-1.1576204) q[3];
sx q[3];
rz(-0.4784085) q[3];
sx q[3];
rz(2.9984409) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
