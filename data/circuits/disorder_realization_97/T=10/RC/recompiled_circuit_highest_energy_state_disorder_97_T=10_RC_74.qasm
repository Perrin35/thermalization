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
rz(-0.47973862) q[0];
sx q[0];
rz(-2.0070183) q[0];
sx q[0];
rz(-1.467508) q[0];
rz(-1.462734) q[1];
sx q[1];
rz(-0.94574133) q[1];
sx q[1];
rz(2.6374964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60668463) q[0];
sx q[0];
rz(-2.1143171) q[0];
sx q[0];
rz(0.96291079) q[0];
rz(-2.7934876) q[2];
sx q[2];
rz(-1.7278302) q[2];
sx q[2];
rz(0.43660313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2254653) q[1];
sx q[1];
rz(-0.91109768) q[1];
sx q[1];
rz(0.74972357) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0023255) q[3];
sx q[3];
rz(-0.33815171) q[3];
sx q[3];
rz(-0.6247181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7291339) q[2];
sx q[2];
rz(-1.2428186) q[2];
sx q[2];
rz(-0.18623713) q[2];
rz(-1.3650182) q[3];
sx q[3];
rz(-0.61187196) q[3];
sx q[3];
rz(1.2046825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81011009) q[0];
sx q[0];
rz(-2.6832566) q[0];
sx q[0];
rz(-0.052074281) q[0];
rz(1.0366108) q[1];
sx q[1];
rz(-0.60984817) q[1];
sx q[1];
rz(-0.61526543) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2589073) q[0];
sx q[0];
rz(-1.8947766) q[0];
sx q[0];
rz(2.2768904) q[0];
rz(-2.6411166) q[2];
sx q[2];
rz(-0.26929528) q[2];
sx q[2];
rz(3.0840906) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76814684) q[1];
sx q[1];
rz(-2.2229338) q[1];
sx q[1];
rz(3.0708205) q[1];
rz(-2.3637487) q[3];
sx q[3];
rz(-1.420971) q[3];
sx q[3];
rz(0.37096393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44876862) q[2];
sx q[2];
rz(-1.8961467) q[2];
sx q[2];
rz(-1.2358865) q[2];
rz(1.1290733) q[3];
sx q[3];
rz(-0.90714199) q[3];
sx q[3];
rz(0.83704078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-1.6579984) q[0];
sx q[0];
rz(-2.3880385) q[0];
sx q[0];
rz(-1.765522) q[0];
rz(-0.54028571) q[1];
sx q[1];
rz(-1.0397725) q[1];
sx q[1];
rz(-1.6859863) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42403635) q[0];
sx q[0];
rz(-0.40923318) q[0];
sx q[0];
rz(-0.17828973) q[0];
rz(-pi) q[1];
rz(2.9380581) q[2];
sx q[2];
rz(-1.4829205) q[2];
sx q[2];
rz(2.3661302) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.27764717) q[1];
sx q[1];
rz(-1.9233812) q[1];
sx q[1];
rz(-0.90971281) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7330879) q[3];
sx q[3];
rz(-2.2271379) q[3];
sx q[3];
rz(1.0711627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49440631) q[2];
sx q[2];
rz(-0.33438412) q[2];
sx q[2];
rz(-1.6953267) q[2];
rz(-1.9485731) q[3];
sx q[3];
rz(-2.1665067) q[3];
sx q[3];
rz(-1.4421991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0874264) q[0];
sx q[0];
rz(-0.26987258) q[0];
sx q[0];
rz(0.42359459) q[0];
rz(2.1845747) q[1];
sx q[1];
rz(-1.7901763) q[1];
sx q[1];
rz(3.0439923) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070484249) q[0];
sx q[0];
rz(-1.1997265) q[0];
sx q[0];
rz(-0.563693) q[0];
rz(-2.9984442) q[2];
sx q[2];
rz(-0.88298702) q[2];
sx q[2];
rz(-1.0873512) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3279539) q[1];
sx q[1];
rz(-2.6662152) q[1];
sx q[1];
rz(2.3085576) q[1];
rz(-pi) q[2];
rz(-0.72928263) q[3];
sx q[3];
rz(-1.3137307) q[3];
sx q[3];
rz(-2.2616539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5264954) q[2];
sx q[2];
rz(-0.47250938) q[2];
sx q[2];
rz(2.9769767) q[2];
rz(-0.047867157) q[3];
sx q[3];
rz(-1.4048301) q[3];
sx q[3];
rz(-2.3544748) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98016244) q[0];
sx q[0];
rz(-0.4466559) q[0];
sx q[0];
rz(1.5572146) q[0];
rz(2.7730675) q[1];
sx q[1];
rz(-1.3216113) q[1];
sx q[1];
rz(-1.6065067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6958862) q[0];
sx q[0];
rz(-2.0793756) q[0];
sx q[0];
rz(2.5832547) q[0];
x q[1];
rz(-2.755411) q[2];
sx q[2];
rz(-2.2791822) q[2];
sx q[2];
rz(2.2798373) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6905744) q[1];
sx q[1];
rz(-1.3214045) q[1];
sx q[1];
rz(1.0728157) q[1];
rz(1.8593349) q[3];
sx q[3];
rz(-1.6464982) q[3];
sx q[3];
rz(2.5385365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2313472) q[2];
sx q[2];
rz(-0.37800899) q[2];
sx q[2];
rz(-2.0965915) q[2];
rz(-2.9735978) q[3];
sx q[3];
rz(-1.955227) q[3];
sx q[3];
rz(-2.7872938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0198233) q[0];
sx q[0];
rz(-1.0291809) q[0];
sx q[0];
rz(-1.2044915) q[0];
rz(2.3380741) q[1];
sx q[1];
rz(-2.3280227) q[1];
sx q[1];
rz(0.71570754) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2432118) q[0];
sx q[0];
rz(-1.0168494) q[0];
sx q[0];
rz(-0.57147632) q[0];
x q[1];
rz(1.7693232) q[2];
sx q[2];
rz(-1.1305222) q[2];
sx q[2];
rz(-2.0581051) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5454272) q[1];
sx q[1];
rz(-1.6837032) q[1];
sx q[1];
rz(-0.66408709) q[1];
x q[2];
rz(-2.0296069) q[3];
sx q[3];
rz(-1.306064) q[3];
sx q[3];
rz(-0.88065016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7314926) q[2];
sx q[2];
rz(-2.4918753) q[2];
sx q[2];
rz(2.0984207) q[2];
rz(0.10797524) q[3];
sx q[3];
rz(-2.2004746) q[3];
sx q[3];
rz(-2.030355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1825948) q[0];
sx q[0];
rz(-1.7549055) q[0];
sx q[0];
rz(1.9166272) q[0];
rz(0.23172465) q[1];
sx q[1];
rz(-2.2544315) q[1];
sx q[1];
rz(-1.8909594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62364336) q[0];
sx q[0];
rz(-0.8276437) q[0];
sx q[0];
rz(-1.736899) q[0];
rz(1.3725946) q[2];
sx q[2];
rz(-1.5491252) q[2];
sx q[2];
rz(-0.54474165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6078203) q[1];
sx q[1];
rz(-1.269286) q[1];
sx q[1];
rz(1.34006) q[1];
rz(-pi) q[2];
rz(1.5579079) q[3];
sx q[3];
rz(-2.1272214) q[3];
sx q[3];
rz(2.9146359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.072824868) q[2];
sx q[2];
rz(-0.67049694) q[2];
sx q[2];
rz(3.0094299) q[2];
rz(-0.69563785) q[3];
sx q[3];
rz(-1.9591103) q[3];
sx q[3];
rz(2.3343991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19145963) q[0];
sx q[0];
rz(-1.5810672) q[0];
sx q[0];
rz(2.4323442) q[0];
rz(-1.6356155) q[1];
sx q[1];
rz(-1.154107) q[1];
sx q[1];
rz(1.0848612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5800672) q[0];
sx q[0];
rz(-0.97082061) q[0];
sx q[0];
rz(-1.6404433) q[0];
x q[1];
rz(0.16895825) q[2];
sx q[2];
rz(-2.3297133) q[2];
sx q[2];
rz(-0.66056992) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0604531) q[1];
sx q[1];
rz(-2.6067002) q[1];
sx q[1];
rz(2.6263531) q[1];
rz(-1.3832757) q[3];
sx q[3];
rz(-0.061358364) q[3];
sx q[3];
rz(-0.96422577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52013493) q[2];
sx q[2];
rz(-1.0026714) q[2];
sx q[2];
rz(-1.3168859) q[2];
rz(0.78553158) q[3];
sx q[3];
rz(-1.1979878) q[3];
sx q[3];
rz(-0.42640105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.921628) q[0];
sx q[0];
rz(-1.4861318) q[0];
sx q[0];
rz(-0.40837902) q[0];
rz(0.95868239) q[1];
sx q[1];
rz(-0.32902333) q[1];
sx q[1];
rz(-1.5214517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085039728) q[0];
sx q[0];
rz(-1.805575) q[0];
sx q[0];
rz(1.1961557) q[0];
x q[1];
rz(1.5217811) q[2];
sx q[2];
rz(-2.1398126) q[2];
sx q[2];
rz(0.70876497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16220763) q[1];
sx q[1];
rz(-1.1623135) q[1];
sx q[1];
rz(0.86345203) q[1];
x q[2];
rz(-1.8835041) q[3];
sx q[3];
rz(-1.7527448) q[3];
sx q[3];
rz(-3.0489717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7659605) q[2];
sx q[2];
rz(-1.1124632) q[2];
sx q[2];
rz(-0.56606236) q[2];
rz(-1.798897) q[3];
sx q[3];
rz(-2.7261966) q[3];
sx q[3];
rz(0.28877637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.600243) q[0];
sx q[0];
rz(-0.039027795) q[0];
sx q[0];
rz(1.716123) q[0];
rz(1.0936945) q[1];
sx q[1];
rz(-0.92233557) q[1];
sx q[1];
rz(-0.38965449) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27586994) q[0];
sx q[0];
rz(-1.8267858) q[0];
sx q[0];
rz(1.1317568) q[0];
x q[1];
rz(-0.25954982) q[2];
sx q[2];
rz(-2.1786111) q[2];
sx q[2];
rz(-1.9454625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13481465) q[1];
sx q[1];
rz(-1.9962427) q[1];
sx q[1];
rz(0.025789217) q[1];
rz(-pi) q[2];
rz(-1.8752619) q[3];
sx q[3];
rz(-1.0926477) q[3];
sx q[3];
rz(-1.4451417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9044372) q[2];
sx q[2];
rz(-0.5390141) q[2];
sx q[2];
rz(-1.944444) q[2];
rz(0.14687471) q[3];
sx q[3];
rz(-0.67320383) q[3];
sx q[3];
rz(-2.1541514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.420153) q[0];
sx q[0];
rz(-2.1099821) q[0];
sx q[0];
rz(-1.4954062) q[0];
rz(-2.738476) q[1];
sx q[1];
rz(-1.8164201) q[1];
sx q[1];
rz(1.4774189) q[1];
rz(-0.18927903) q[2];
sx q[2];
rz(-0.58313417) q[2];
sx q[2];
rz(-1.492576) q[2];
rz(2.7561989) q[3];
sx q[3];
rz(-1.7447532) q[3];
sx q[3];
rz(0.90739653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
