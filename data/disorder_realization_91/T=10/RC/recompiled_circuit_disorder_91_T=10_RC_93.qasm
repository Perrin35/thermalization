OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4361753) q[0];
sx q[0];
rz(-0.56646148) q[0];
sx q[0];
rz(0.17106549) q[0];
rz(-1.6879727) q[1];
sx q[1];
rz(3.4847335) q[1];
sx q[1];
rz(7.6141678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91141191) q[0];
sx q[0];
rz(-2.4898306) q[0];
sx q[0];
rz(3.0551811) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5924442) q[2];
sx q[2];
rz(-2.0865555) q[2];
sx q[2];
rz(1.6025008) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1316489) q[1];
sx q[1];
rz(-2.8899) q[1];
sx q[1];
rz(-1.5021055) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64294502) q[3];
sx q[3];
rz(-1.7503947) q[3];
sx q[3];
rz(-2.5760866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77582899) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(-1.8120871) q[2];
rz(1.4154411) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(-0.20234385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.99130327) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(-1.6888899) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-1.0849489) q[1];
sx q[1];
rz(0.30028775) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4703003) q[0];
sx q[0];
rz(-1.4503345) q[0];
sx q[0];
rz(-2.5380773) q[0];
x q[1];
rz(-0.57057256) q[2];
sx q[2];
rz(-1.3010446) q[2];
sx q[2];
rz(-2.3938993) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.049388) q[1];
sx q[1];
rz(-0.37282473) q[1];
sx q[1];
rz(-0.91918175) q[1];
x q[2];
rz(-0.16365666) q[3];
sx q[3];
rz(-1.7573342) q[3];
sx q[3];
rz(-1.2638307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.369027) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(1.0380113) q[2];
rz(-2.299262) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(2.7712908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5623986) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(-2.753479) q[0];
rz(-0.072470486) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(0.31633502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65773327) q[0];
sx q[0];
rz(-1.8158001) q[0];
sx q[0];
rz(0.0018951456) q[0];
x q[1];
rz(-0.57245589) q[2];
sx q[2];
rz(-1.2081895) q[2];
sx q[2];
rz(0.87871694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9005147) q[1];
sx q[1];
rz(-1.1763651) q[1];
sx q[1];
rz(-2.6532252) q[1];
x q[2];
rz(1.9100788) q[3];
sx q[3];
rz(-1.2959891) q[3];
sx q[3];
rz(0.85698444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0324273) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(0.16080984) q[2];
rz(3.0155449) q[3];
sx q[3];
rz(-1.3879317) q[3];
sx q[3];
rz(1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9008824) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(2.3098992) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(-1.5136738) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31878372) q[0];
sx q[0];
rz(-1.6167287) q[0];
sx q[0];
rz(2.4687693) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1349147) q[2];
sx q[2];
rz(-1.5022105) q[2];
sx q[2];
rz(0.90142957) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6935389) q[1];
sx q[1];
rz(-1.5539907) q[1];
sx q[1];
rz(2.0190713) q[1];
x q[2];
rz(1.1553331) q[3];
sx q[3];
rz(-0.88298015) q[3];
sx q[3];
rz(-0.26879877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4251129) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(-1.9173737) q[2];
rz(-0.41695693) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.083374627) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(1.303724) q[0];
rz(2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(0.051503332) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542959) q[0];
sx q[0];
rz(-1.686839) q[0];
sx q[0];
rz(-0.10033484) q[0];
x q[1];
rz(-2.0578458) q[2];
sx q[2];
rz(-0.32099989) q[2];
sx q[2];
rz(-3.0082862) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3225973) q[1];
sx q[1];
rz(-1.4154735) q[1];
sx q[1];
rz(1.1197234) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0052161) q[3];
sx q[3];
rz(-2.4032421) q[3];
sx q[3];
rz(-2.7969558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45433989) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(0.59801897) q[2];
rz(-0.55001843) q[3];
sx q[3];
rz(-0.23967448) q[3];
sx q[3];
rz(-1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0710058) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(-2.9396074) q[0];
rz(0.96549353) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(0.083267033) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9986388) q[0];
sx q[0];
rz(-0.26252258) q[0];
sx q[0];
rz(0.72857626) q[0];
rz(-pi) q[1];
rz(-2.3699049) q[2];
sx q[2];
rz(-0.37479127) q[2];
sx q[2];
rz(0.37116606) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0373842) q[1];
sx q[1];
rz(-1.6777778) q[1];
sx q[1];
rz(-1.2030829) q[1];
rz(2.5842651) q[3];
sx q[3];
rz(-1.8263655) q[3];
sx q[3];
rz(2.717358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0319556) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(-1.7588245) q[2];
rz(-0.2494732) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270585) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(-0.89685857) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-2.9607594) q[1];
sx q[1];
rz(2.0239963) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5516978) q[0];
sx q[0];
rz(-2.4015744) q[0];
sx q[0];
rz(0.19703534) q[0];
rz(-pi) q[1];
rz(-2.2824077) q[2];
sx q[2];
rz(-1.7297598) q[2];
sx q[2];
rz(-0.32260103) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7270131) q[1];
sx q[1];
rz(-1.0158744) q[1];
sx q[1];
rz(-1.8068061) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8459122) q[3];
sx q[3];
rz(-1.0529622) q[3];
sx q[3];
rz(0.3412316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0573132) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(-0.21952595) q[2];
rz(-2.7548742) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893843) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(2.9242933) q[0];
rz(-3.1106588) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(-0.6589748) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17651672) q[0];
sx q[0];
rz(-2.342431) q[0];
sx q[0];
rz(-1.6060711) q[0];
rz(-pi) q[1];
rz(0.89491567) q[2];
sx q[2];
rz(-0.96940982) q[2];
sx q[2];
rz(-0.88047699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1606071) q[1];
sx q[1];
rz(-1.7874582) q[1];
sx q[1];
rz(-2.3807081) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6625255) q[3];
sx q[3];
rz(-0.64507161) q[3];
sx q[3];
rz(-2.0446387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.7610901) q[2];
sx q[2];
rz(2.8213815) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.1956918) q[3];
sx q[3];
rz(1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76072389) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(-2.4587801) q[0];
rz(2.071351) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(1.210093) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4142128) q[0];
sx q[0];
rz(-2.7276037) q[0];
sx q[0];
rz(-1.0340286) q[0];
rz(-pi) q[1];
rz(2.1754873) q[2];
sx q[2];
rz(-1.7793057) q[2];
sx q[2];
rz(1.0143806) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2883816) q[1];
sx q[1];
rz(-2.2182811) q[1];
sx q[1];
rz(-2.7700469) q[1];
rz(-1.1990511) q[3];
sx q[3];
rz(-0.98402714) q[3];
sx q[3];
rz(1.1128807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5658297) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(-2.5749717) q[2];
rz(0.9295272) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(-1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(1.7096827) q[0];
rz(-2.6152949) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(0.79968232) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4096217) q[0];
sx q[0];
rz(-0.35847607) q[0];
sx q[0];
rz(1.8147857) q[0];
rz(0.16074796) q[2];
sx q[2];
rz(-1.4904009) q[2];
sx q[2];
rz(2.423167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29337063) q[1];
sx q[1];
rz(-1.231041) q[1];
sx q[1];
rz(2.1128224) q[1];
x q[2];
rz(1.6465854) q[3];
sx q[3];
rz(-1.4255376) q[3];
sx q[3];
rz(2.4491572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8563103) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(2.4463859) q[2];
rz(2.5837512) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.2790537) q[1];
sx q[1];
rz(-0.74395724) q[1];
sx q[1];
rz(-0.67768135) q[1];
rz(1.5699301) q[2];
sx q[2];
rz(-1.2348839) q[2];
sx q[2];
rz(0.15765794) q[2];
rz(2.9813319) q[3];
sx q[3];
rz(-2.3370623) q[3];
sx q[3];
rz(-0.2094895) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
