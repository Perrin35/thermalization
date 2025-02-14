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
rz(0.12240527) q[0];
sx q[0];
rz(2.2221017) q[0];
sx q[0];
rz(10.623951) q[0];
rz(0.1872669) q[1];
sx q[1];
rz(-2.5993102) q[1];
sx q[1];
rz(-1.5195001) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6698201) q[0];
sx q[0];
rz(-2.5216521) q[0];
sx q[0];
rz(0.13373904) q[0];
rz(-pi) q[1];
rz(-0.3737195) q[2];
sx q[2];
rz(-1.1454095) q[2];
sx q[2];
rz(-0.56217867) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70473814) q[1];
sx q[1];
rz(-0.66009843) q[1];
sx q[1];
rz(2.4515602) q[1];
rz(1.5222129) q[3];
sx q[3];
rz(-1.8147) q[3];
sx q[3];
rz(2.3590306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2970994) q[2];
sx q[2];
rz(-2.5145734) q[2];
sx q[2];
rz(-1.373488) q[2];
rz(-2.3376236) q[3];
sx q[3];
rz(-1.6425902) q[3];
sx q[3];
rz(1.1871626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3961769) q[0];
sx q[0];
rz(-2.4276623) q[0];
sx q[0];
rz(-0.3748689) q[0];
rz(-0.51775852) q[1];
sx q[1];
rz(-1.9381783) q[1];
sx q[1];
rz(1.7727324) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8086209) q[0];
sx q[0];
rz(-3.140641) q[0];
sx q[0];
rz(-0.10097058) q[0];
rz(-pi) q[1];
rz(-2.9451319) q[2];
sx q[2];
rz(-1.5443373) q[2];
sx q[2];
rz(-0.65633869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9437286) q[1];
sx q[1];
rz(-1.1625338) q[1];
sx q[1];
rz(0.54710435) q[1];
rz(-3.0053455) q[3];
sx q[3];
rz(-1.7833424) q[3];
sx q[3];
rz(3.1177136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0181197) q[2];
sx q[2];
rz(-0.70088434) q[2];
sx q[2];
rz(2.4269721) q[2];
rz(2.9361652) q[3];
sx q[3];
rz(-1.8193865) q[3];
sx q[3];
rz(-1.2986758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4648723) q[0];
sx q[0];
rz(-2.1045852) q[0];
sx q[0];
rz(0.49764693) q[0];
rz(-0.63703713) q[1];
sx q[1];
rz(-0.62890816) q[1];
sx q[1];
rz(3.0348437) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9861988) q[0];
sx q[0];
rz(-1.2717016) q[0];
sx q[0];
rz(0.22217447) q[0];
rz(-1.5791513) q[2];
sx q[2];
rz(-0.66973493) q[2];
sx q[2];
rz(1.005203) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6124668) q[1];
sx q[1];
rz(-0.71536109) q[1];
sx q[1];
rz(1.4892519) q[1];
x q[2];
rz(2.51744) q[3];
sx q[3];
rz(-0.48138967) q[3];
sx q[3];
rz(2.2311051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5367624) q[2];
sx q[2];
rz(-2.3894775) q[2];
sx q[2];
rz(2.9467648) q[2];
rz(-0.005006494) q[3];
sx q[3];
rz(-0.61401335) q[3];
sx q[3];
rz(-2.3495638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.1057338) q[0];
sx q[0];
rz(-0.53426131) q[0];
sx q[0];
rz(2.1015097) q[0];
rz(0.19373521) q[1];
sx q[1];
rz(-1.5015142) q[1];
sx q[1];
rz(0.74660444) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6116007) q[0];
sx q[0];
rz(-1.1318637) q[0];
sx q[0];
rz(-0.12719391) q[0];
x q[1];
rz(0.87073054) q[2];
sx q[2];
rz(-1.680003) q[2];
sx q[2];
rz(2.2094215) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6317217) q[1];
sx q[1];
rz(-1.5190017) q[1];
sx q[1];
rz(-2.447261) q[1];
x q[2];
rz(1.6799404) q[3];
sx q[3];
rz(-1.3323083) q[3];
sx q[3];
rz(1.3772688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1497583) q[2];
sx q[2];
rz(-1.4703625) q[2];
sx q[2];
rz(0.20450083) q[2];
rz(1.9483942) q[3];
sx q[3];
rz(-2.2843993) q[3];
sx q[3];
rz(-2.1804555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0123154) q[0];
sx q[0];
rz(-0.17450541) q[0];
sx q[0];
rz(-1.0804863) q[0];
rz(-0.37733817) q[1];
sx q[1];
rz(-1.6308547) q[1];
sx q[1];
rz(-0.89471716) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52496929) q[0];
sx q[0];
rz(-1.3436514) q[0];
sx q[0];
rz(0.36360111) q[0];
x q[1];
rz(0.65308833) q[2];
sx q[2];
rz(-0.90735596) q[2];
sx q[2];
rz(1.1764248) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.0061568312) q[1];
sx q[1];
rz(-0.73169152) q[1];
sx q[1];
rz(-2.3680971) q[1];
rz(0.36578806) q[3];
sx q[3];
rz(-2.180763) q[3];
sx q[3];
rz(-0.47780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6696024) q[2];
sx q[2];
rz(-2.695638) q[2];
sx q[2];
rz(-0.10776821) q[2];
rz(-2.9660411) q[3];
sx q[3];
rz(-1.8864417) q[3];
sx q[3];
rz(0.56041437) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948697) q[0];
sx q[0];
rz(-0.51114285) q[0];
sx q[0];
rz(1.0119337) q[0];
rz(1.5049505) q[1];
sx q[1];
rz(-2.4220146) q[1];
sx q[1];
rz(-1.0171657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95421529) q[0];
sx q[0];
rz(-0.34437505) q[0];
sx q[0];
rz(-1.0402753) q[0];
x q[1];
rz(2.1332729) q[2];
sx q[2];
rz(-0.89300821) q[2];
sx q[2];
rz(-0.00046367292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.519327) q[1];
sx q[1];
rz(-0.19058386) q[1];
sx q[1];
rz(-3.0960073) q[1];
rz(1.1215707) q[3];
sx q[3];
rz(-1.7127081) q[3];
sx q[3];
rz(2.6117539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4096058) q[2];
sx q[2];
rz(-0.66183949) q[2];
sx q[2];
rz(-0.99676639) q[2];
rz(0.064920001) q[3];
sx q[3];
rz(-1.7305948) q[3];
sx q[3];
rz(1.9916649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0885334) q[0];
sx q[0];
rz(-0.94045883) q[0];
sx q[0];
rz(0.9084107) q[0];
rz(2.9216595) q[1];
sx q[1];
rz(-1.5989774) q[1];
sx q[1];
rz(-1.251108) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.096643) q[0];
sx q[0];
rz(-1.310955) q[0];
sx q[0];
rz(-2.0578865) q[0];
rz(-pi) q[1];
x q[1];
rz(1.471453) q[2];
sx q[2];
rz(-1.8548428) q[2];
sx q[2];
rz(2.5104475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0617705) q[1];
sx q[1];
rz(-2.1822189) q[1];
sx q[1];
rz(-2.3794214) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7102107) q[3];
sx q[3];
rz(-2.406139) q[3];
sx q[3];
rz(-2.9289587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1088341) q[2];
sx q[2];
rz(-2.3075576) q[2];
sx q[2];
rz(-1.533482) q[2];
rz(-1.1533302) q[3];
sx q[3];
rz(-1.0554375) q[3];
sx q[3];
rz(-1.4920894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9913427) q[0];
sx q[0];
rz(-0.38337502) q[0];
sx q[0];
rz(0.94069329) q[0];
rz(-0.13537814) q[1];
sx q[1];
rz(-1.7372513) q[1];
sx q[1];
rz(1.0248331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.301046) q[0];
sx q[0];
rz(-2.4468166) q[0];
sx q[0];
rz(-1.077681) q[0];
rz(1.631279) q[2];
sx q[2];
rz(-0.96394682) q[2];
sx q[2];
rz(-1.902193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.35541818) q[1];
sx q[1];
rz(-2.9045331) q[1];
sx q[1];
rz(2.5393344) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5190184) q[3];
sx q[3];
rz(-1.9815784) q[3];
sx q[3];
rz(0.67922986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3465053) q[2];
sx q[2];
rz(-2.2515209) q[2];
sx q[2];
rz(-0.65004641) q[2];
rz(2.5551689) q[3];
sx q[3];
rz(-1.6610049) q[3];
sx q[3];
rz(0.75064269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8745678) q[0];
sx q[0];
rz(-2.6331007) q[0];
sx q[0];
rz(-1.9961927) q[0];
rz(1.3151431) q[1];
sx q[1];
rz(-1.2329085) q[1];
sx q[1];
rz(2.4536536) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8696339) q[0];
sx q[0];
rz(-2.7209073) q[0];
sx q[0];
rz(-1.0933594) q[0];
x q[1];
rz(-0.97369377) q[2];
sx q[2];
rz(-1.2427689) q[2];
sx q[2];
rz(-0.51842481) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6932144) q[1];
sx q[1];
rz(-0.82469392) q[1];
sx q[1];
rz(-1.34971) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93291544) q[3];
sx q[3];
rz(-2.5543) q[3];
sx q[3];
rz(1.4838532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73100662) q[2];
sx q[2];
rz(-1.1062016) q[2];
sx q[2];
rz(-2.2753687) q[2];
rz(0.90803641) q[3];
sx q[3];
rz(-2.3767545) q[3];
sx q[3];
rz(-1.0959371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1753801) q[0];
sx q[0];
rz(-1.9575653) q[0];
sx q[0];
rz(-0.41561919) q[0];
rz(0.79187727) q[1];
sx q[1];
rz(-1.917058) q[1];
sx q[1];
rz(-0.22722879) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4675605) q[0];
sx q[0];
rz(-1.5825543) q[0];
sx q[0];
rz(-3.1255162) q[0];
rz(-0.80663075) q[2];
sx q[2];
rz(-1.5147594) q[2];
sx q[2];
rz(-3.0174668) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91252226) q[1];
sx q[1];
rz(-2.4991112) q[1];
sx q[1];
rz(0.43825348) q[1];
rz(-pi) q[2];
rz(0.51550224) q[3];
sx q[3];
rz(-1.1519264) q[3];
sx q[3];
rz(-0.91332289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26433429) q[2];
sx q[2];
rz(-2.4206968) q[2];
sx q[2];
rz(-0.43992511) q[2];
rz(2.544493) q[3];
sx q[3];
rz(-0.66949451) q[3];
sx q[3];
rz(1.7841024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6482342) q[0];
sx q[0];
rz(-1.5536722) q[0];
sx q[0];
rz(-1.8552725) q[0];
rz(-1.413912) q[1];
sx q[1];
rz(-1.021011) q[1];
sx q[1];
rz(0.11722142) q[1];
rz(1.5154503) q[2];
sx q[2];
rz(-1.9774441) q[2];
sx q[2];
rz(-1.0033506) q[2];
rz(0.24735484) q[3];
sx q[3];
rz(-0.83360278) q[3];
sx q[3];
rz(-3.0616888) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
