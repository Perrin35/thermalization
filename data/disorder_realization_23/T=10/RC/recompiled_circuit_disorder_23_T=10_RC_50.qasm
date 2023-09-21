OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8795348) q[0];
sx q[0];
rz(-1.4095925) q[0];
sx q[0];
rz(-1.4341266) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.943676) q[0];
sx q[0];
rz(-1.1063873) q[0];
sx q[0];
rz(0.29505131) q[0];
rz(-pi) q[1];
rz(-0.11169545) q[2];
sx q[2];
rz(-1.4354435) q[2];
sx q[2];
rz(0.40458194) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5814372) q[1];
sx q[1];
rz(-1.1168224) q[1];
sx q[1];
rz(-0.27172383) q[1];
rz(-pi) q[2];
rz(2.0166287) q[3];
sx q[3];
rz(-2.6615407) q[3];
sx q[3];
rz(-1.1207086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3216386) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(-2.3036172) q[2];
rz(-0.49301246) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943966) q[0];
sx q[0];
rz(-0.72421873) q[0];
sx q[0];
rz(-1.2778506) q[0];
rz(-0.17678075) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(2.7094254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1457739) q[0];
sx q[0];
rz(-2.5403025) q[0];
sx q[0];
rz(2.6390618) q[0];
rz(-pi) q[1];
rz(1.997666) q[2];
sx q[2];
rz(-1.9383213) q[2];
sx q[2];
rz(1.2183684) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9464311) q[1];
sx q[1];
rz(-2.6841607) q[1];
sx q[1];
rz(-1.3034348) q[1];
rz(-pi) q[2];
rz(2.1357972) q[3];
sx q[3];
rz(-1.5917935) q[3];
sx q[3];
rz(2.3137623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54923487) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(-0.48669997) q[2];
rz(-1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(-0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.10087) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(2.7242463) q[0];
rz(1.4886645) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(0.506385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15784141) q[0];
sx q[0];
rz(-2.7451773) q[0];
sx q[0];
rz(-2.0435964) q[0];
rz(-0.26489139) q[2];
sx q[2];
rz(-0.75220097) q[2];
sx q[2];
rz(2.0843992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8295146) q[1];
sx q[1];
rz(-1.8784338) q[1];
sx q[1];
rz(2.679146) q[1];
rz(-pi) q[2];
rz(-2.9590928) q[3];
sx q[3];
rz(-1.2068519) q[3];
sx q[3];
rz(2.5516627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4425519) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(0.59147269) q[2];
rz(-0.58602035) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1699003) q[0];
sx q[0];
rz(-0.62830347) q[0];
sx q[0];
rz(-0.83918321) q[0];
rz(0.025578586) q[1];
sx q[1];
rz(-2.4459116) q[1];
sx q[1];
rz(1.5930088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4913113) q[0];
sx q[0];
rz(-2.6822753) q[0];
sx q[0];
rz(1.7772872) q[0];
rz(-2.3392802) q[2];
sx q[2];
rz(-2.1244086) q[2];
sx q[2];
rz(0.72788903) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5316481) q[1];
sx q[1];
rz(-0.87024401) q[1];
sx q[1];
rz(-1.0611666) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43058331) q[3];
sx q[3];
rz(-1.8292556) q[3];
sx q[3];
rz(-2.5653198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8362391) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(-0.099686064) q[2];
rz(-2.1827407) q[3];
sx q[3];
rz(-1.3189664) q[3];
sx q[3];
rz(-1.3249741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56617671) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(-0.25594041) q[0];
rz(0.4610962) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(2.3815313) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5998659) q[0];
sx q[0];
rz(-1.6875629) q[0];
sx q[0];
rz(-1.6745425) q[0];
rz(1.9854457) q[2];
sx q[2];
rz(-1.6290602) q[2];
sx q[2];
rz(2.0410048) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0415503) q[1];
sx q[1];
rz(-1.4428382) q[1];
sx q[1];
rz(-3.0138739) q[1];
rz(-pi) q[2];
rz(-1.101196) q[3];
sx q[3];
rz(-2.4393775) q[3];
sx q[3];
rz(1.5839674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5710859) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(0.6742397) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-2.6847697) q[3];
sx q[3];
rz(3.1242483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6102585) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(-1.9624788) q[0];
rz(0.20482652) q[1];
sx q[1];
rz(-0.79524672) q[1];
sx q[1];
rz(-2.0746322) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.897192) q[0];
sx q[0];
rz(-3.1163437) q[0];
sx q[0];
rz(-0.61141725) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1322137) q[2];
sx q[2];
rz(-1.9595993) q[2];
sx q[2];
rz(-2.4593381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3031591) q[1];
sx q[1];
rz(-1.5232956) q[1];
sx q[1];
rz(-0.82364239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9818929) q[3];
sx q[3];
rz(-2.3387863) q[3];
sx q[3];
rz(-1.476895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.431488) q[2];
sx q[2];
rz(-1.2892712) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(-2.4152749) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(-2.8360951) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85754919) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(1.42111) q[0];
rz(-2.9395318) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(0.85817671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7875123) q[0];
sx q[0];
rz(-1.6991827) q[0];
sx q[0];
rz(-1.7187198) q[0];
rz(-1.981295) q[2];
sx q[2];
rz(-1.5294642) q[2];
sx q[2];
rz(-1.873204) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.757829) q[1];
sx q[1];
rz(-1.5997636) q[1];
sx q[1];
rz(1.5555192) q[1];
rz(-pi) q[2];
rz(-1.2248366) q[3];
sx q[3];
rz(-2.1145027) q[3];
sx q[3];
rz(-0.35017761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8537366) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(1.8161592) q[2];
rz(2.251513) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(-0.98852283) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778075) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(-2.7668787) q[0];
rz(-2.162714) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(-1.3495061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50586787) q[0];
sx q[0];
rz(-0.74698193) q[0];
sx q[0];
rz(2.5542459) q[0];
rz(3.1153203) q[2];
sx q[2];
rz(-2.4233344) q[2];
sx q[2];
rz(0.19933137) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6264682) q[1];
sx q[1];
rz(-1.5212458) q[1];
sx q[1];
rz(1.7527761) q[1];
rz(-pi) q[2];
rz(-0.41534822) q[3];
sx q[3];
rz(-0.58598622) q[3];
sx q[3];
rz(-1.8893482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4902041) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(2.6728969) q[2];
rz(-1.9474585) q[3];
sx q[3];
rz(-1.8374551) q[3];
sx q[3];
rz(1.7780001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63012183) q[0];
sx q[0];
rz(-2.2352495) q[0];
sx q[0];
rz(1.8796896) q[0];
rz(-2.966554) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(1.6040241) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14311929) q[0];
sx q[0];
rz(-2.7612918) q[0];
sx q[0];
rz(-0.93912504) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41861694) q[2];
sx q[2];
rz(-1.6659684) q[2];
sx q[2];
rz(-1.0375432) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21737145) q[1];
sx q[1];
rz(-1.4573759) q[1];
sx q[1];
rz(2.4356615) q[1];
x q[2];
rz(-2.2213307) q[3];
sx q[3];
rz(-0.55763054) q[3];
sx q[3];
rz(-2.695431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.039915446) q[2];
sx q[2];
rz(-0.95255178) q[2];
sx q[2];
rz(-0.76134479) q[2];
rz(-2.2411761) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(0.049023978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3392357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(-0.21690579) q[0];
rz(-2.5096109) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(-2.1868618) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2749598) q[0];
sx q[0];
rz(-1.2506335) q[0];
sx q[0];
rz(-2.3812256) q[0];
rz(1.9913313) q[2];
sx q[2];
rz(-1.0610126) q[2];
sx q[2];
rz(1.860294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.875647) q[1];
sx q[1];
rz(-1.6469643) q[1];
sx q[1];
rz(0.45653685) q[1];
x q[2];
rz(-2.6984152) q[3];
sx q[3];
rz(-2.157353) q[3];
sx q[3];
rz(3.0122258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0163991) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(2.1968502) q[2];
rz(0.38481209) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(2.184536) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64086296) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-2.2568933) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(0.83512886) q[2];
sx q[2];
rz(-2.0616812) q[2];
sx q[2];
rz(2.4731935) q[2];
rz(-1.1273884) q[3];
sx q[3];
rz(-1.3848806) q[3];
sx q[3];
rz(1.0564907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
