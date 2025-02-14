OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0394734) q[0];
sx q[0];
rz(-1.4730299) q[0];
sx q[0];
rz(0.13134512) q[0];
rz(-3.1047473) q[1];
sx q[1];
rz(-2.614202) q[1];
sx q[1];
rz(-1.8324469) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1037169) q[0];
sx q[0];
rz(-1.193422) q[0];
sx q[0];
rz(-0.81822936) q[0];
x q[1];
rz(-0.45707656) q[2];
sx q[2];
rz(-2.0684048) q[2];
sx q[2];
rz(3.0534985) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75571698) q[1];
sx q[1];
rz(-0.39744032) q[1];
sx q[1];
rz(0.65544513) q[1];
rz(0.18351002) q[3];
sx q[3];
rz(-1.7391296) q[3];
sx q[3];
rz(-2.3678697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5347791) q[2];
sx q[2];
rz(-0.35627347) q[2];
sx q[2];
rz(0.66205364) q[2];
rz(-0.74696294) q[3];
sx q[3];
rz(-2.2253939) q[3];
sx q[3];
rz(-1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58716431) q[0];
sx q[0];
rz(-2.9968408) q[0];
sx q[0];
rz(-1.9594877) q[0];
rz(2.3960522) q[1];
sx q[1];
rz(-0.59919557) q[1];
sx q[1];
rz(0.10339698) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.11907) q[0];
sx q[0];
rz(-1.8028717) q[0];
sx q[0];
rz(-2.705535) q[0];
rz(-0.38455268) q[2];
sx q[2];
rz(-1.0083303) q[2];
sx q[2];
rz(2.1484274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2491746) q[1];
sx q[1];
rz(-1.2399481) q[1];
sx q[1];
rz(0.82280092) q[1];
rz(-2.313226) q[3];
sx q[3];
rz(-1.4423372) q[3];
sx q[3];
rz(2.1118233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5891002) q[2];
sx q[2];
rz(-1.2673667) q[2];
sx q[2];
rz(2.6972771) q[2];
rz(2.0878504) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(-1.1963371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0788197) q[0];
sx q[0];
rz(-1.2508996) q[0];
sx q[0];
rz(2.479082) q[0];
rz(-1.6569116) q[1];
sx q[1];
rz(-1.0228446) q[1];
sx q[1];
rz(0.2167162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0652567) q[0];
sx q[0];
rz(-2.2562593) q[0];
sx q[0];
rz(0.61409593) q[0];
x q[1];
rz(-2.6609801) q[2];
sx q[2];
rz(-0.79250249) q[2];
sx q[2];
rz(-2.8414244) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.81074508) q[1];
sx q[1];
rz(-2.243089) q[1];
sx q[1];
rz(-0.051203392) q[1];
x q[2];
rz(-1.4101348) q[3];
sx q[3];
rz(-1.0815797) q[3];
sx q[3];
rz(2.1833724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4776939) q[2];
sx q[2];
rz(-2.615216) q[2];
sx q[2];
rz(-2.9065175) q[2];
rz(2.7663686) q[3];
sx q[3];
rz(-2.2347361) q[3];
sx q[3];
rz(2.4231329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5594056) q[0];
sx q[0];
rz(-3.0538054) q[0];
sx q[0];
rz(1.0002332) q[0];
rz(1.2031215) q[1];
sx q[1];
rz(-1.6327881) q[1];
sx q[1];
rz(2.5968754) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9207912) q[0];
sx q[0];
rz(-1.5259597) q[0];
sx q[0];
rz(2.4283571) q[0];
x q[1];
rz(-2.4248872) q[2];
sx q[2];
rz(-2.1726328) q[2];
sx q[2];
rz(-0.81317893) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1511606) q[1];
sx q[1];
rz(-2.9176788) q[1];
sx q[1];
rz(-2.1316705) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9849378) q[3];
sx q[3];
rz(-2.4865287) q[3];
sx q[3];
rz(-2.6998346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7228221) q[2];
sx q[2];
rz(-2.3493769) q[2];
sx q[2];
rz(2.344632) q[2];
rz(0.289251) q[3];
sx q[3];
rz(-2.1713493) q[3];
sx q[3];
rz(2.7204035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61409426) q[0];
sx q[0];
rz(-0.088936381) q[0];
sx q[0];
rz(1.2689137) q[0];
rz(-2.1753963) q[1];
sx q[1];
rz(-1.7485917) q[1];
sx q[1];
rz(2.6752245) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.120935) q[0];
sx q[0];
rz(-1.0758721) q[0];
sx q[0];
rz(2.0536971) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5086002) q[2];
sx q[2];
rz(-2.5883753) q[2];
sx q[2];
rz(-2.4621682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1524017) q[1];
sx q[1];
rz(-0.70211239) q[1];
sx q[1];
rz(1.0662119) q[1];
rz(0.90762805) q[3];
sx q[3];
rz(-1.25718) q[3];
sx q[3];
rz(0.37330353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7352778) q[2];
sx q[2];
rz(-0.80253989) q[2];
sx q[2];
rz(-0.80424133) q[2];
rz(-2.6546226) q[3];
sx q[3];
rz(-2.099497) q[3];
sx q[3];
rz(0.0074726661) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2020579) q[0];
sx q[0];
rz(-2.845293) q[0];
sx q[0];
rz(-0.95788389) q[0];
rz(-1.6288039) q[1];
sx q[1];
rz(-1.0572546) q[1];
sx q[1];
rz(1.5527976) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88547546) q[0];
sx q[0];
rz(-1.6506356) q[0];
sx q[0];
rz(-3.1186597) q[0];
rz(-0.79926305) q[2];
sx q[2];
rz(-2.7529319) q[2];
sx q[2];
rz(2.6262019) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.54974906) q[1];
sx q[1];
rz(-1.3938245) q[1];
sx q[1];
rz(0.66284499) q[1];
rz(-pi) q[2];
rz(1.1670052) q[3];
sx q[3];
rz(-0.98282114) q[3];
sx q[3];
rz(2.605632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2019041) q[2];
sx q[2];
rz(-2.0857911) q[2];
sx q[2];
rz(-2.4064257) q[2];
rz(0.9489263) q[3];
sx q[3];
rz(-2.2831423) q[3];
sx q[3];
rz(0.86400509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76148024) q[0];
sx q[0];
rz(-1.004847) q[0];
sx q[0];
rz(0.6066221) q[0];
rz(2.3573549) q[1];
sx q[1];
rz(-1.3332858) q[1];
sx q[1];
rz(1.7971136) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8356268) q[0];
sx q[0];
rz(-2.2971662) q[0];
sx q[0];
rz(2.4854922) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8227578) q[2];
sx q[2];
rz(-2.0540385) q[2];
sx q[2];
rz(-1.9117219) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.575702) q[1];
sx q[1];
rz(-1.9676349) q[1];
sx q[1];
rz(1.3299408) q[1];
rz(-pi) q[2];
rz(-1.3913888) q[3];
sx q[3];
rz(-0.47166892) q[3];
sx q[3];
rz(0.52667945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.532423) q[2];
sx q[2];
rz(-0.50749856) q[2];
sx q[2];
rz(1.3896821) q[2];
rz(2.9621647) q[3];
sx q[3];
rz(-0.98589412) q[3];
sx q[3];
rz(0.40670407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0824025) q[0];
sx q[0];
rz(-2.817509) q[0];
sx q[0];
rz(-0.75505906) q[0];
rz(0.45626196) q[1];
sx q[1];
rz(-1.4249964) q[1];
sx q[1];
rz(1.0770575) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74769831) q[0];
sx q[0];
rz(-2.1222881) q[0];
sx q[0];
rz(1.3026139) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78387129) q[2];
sx q[2];
rz(-1.4848564) q[2];
sx q[2];
rz(2.849907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7491319) q[1];
sx q[1];
rz(-1.8549671) q[1];
sx q[1];
rz(-2.6662213) q[1];
x q[2];
rz(-1.0302587) q[3];
sx q[3];
rz(-1.4344627) q[3];
sx q[3];
rz(0.34878525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72851744) q[2];
sx q[2];
rz(-2.04144) q[2];
sx q[2];
rz(2.4033578) q[2];
rz(1.632656) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(1.9807321) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64410925) q[0];
sx q[0];
rz(-1.6832385) q[0];
sx q[0];
rz(2.4926376) q[0];
rz(-2.451918) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(-2.7241657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8807424) q[0];
sx q[0];
rz(-2.1584903) q[0];
sx q[0];
rz(-1.8729314) q[0];
rz(-0.242495) q[2];
sx q[2];
rz(-2.8333377) q[2];
sx q[2];
rz(2.2876842) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74655731) q[1];
sx q[1];
rz(-0.42320028) q[1];
sx q[1];
rz(2.8854135) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40006684) q[3];
sx q[3];
rz(-1.2037945) q[3];
sx q[3];
rz(2.7719967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6680341) q[2];
sx q[2];
rz(-1.4213976) q[2];
sx q[2];
rz(-3.0958946) q[2];
rz(-2.8081196) q[3];
sx q[3];
rz(-2.5662751) q[3];
sx q[3];
rz(-0.12588178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68542737) q[0];
sx q[0];
rz(-2.6021155) q[0];
sx q[0];
rz(-1.217655) q[0];
rz(0.24670163) q[1];
sx q[1];
rz(-1.8496937) q[1];
sx q[1];
rz(0.47952476) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96146255) q[0];
sx q[0];
rz(-0.96213522) q[0];
sx q[0];
rz(0.5345688) q[0];
rz(-pi) q[1];
rz(-2.2662451) q[2];
sx q[2];
rz(-1.9672868) q[2];
sx q[2];
rz(0.04405313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.785276) q[1];
sx q[1];
rz(-2.0802167) q[1];
sx q[1];
rz(1.2356349) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8080975) q[3];
sx q[3];
rz(-1.1166683) q[3];
sx q[3];
rz(0.25179201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5852927) q[2];
sx q[2];
rz(-1.8712529) q[2];
sx q[2];
rz(-1.4053819) q[2];
rz(-2.4893238) q[3];
sx q[3];
rz(-2.8758719) q[3];
sx q[3];
rz(-2.4238267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.9115059) q[0];
sx q[0];
rz(-1.390504) q[0];
sx q[0];
rz(2.6934296) q[0];
rz(-2.3975092) q[1];
sx q[1];
rz(-2.4241445) q[1];
sx q[1];
rz(-1.4345899) q[1];
rz(-2.9058331) q[2];
sx q[2];
rz(-2.5452062) q[2];
sx q[2];
rz(1.7231981) q[2];
rz(-0.99540972) q[3];
sx q[3];
rz(-1.8707471) q[3];
sx q[3];
rz(-1.7729014) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
